# 真實基線為 exp：生成、配適與動態 Brier Score

研究設計為生成一份真實基線為 exponential 的資料，再使用同一份資料配適不同聯合模型，比較固定預測跨度下的動態 Brier score。生成端採用有限追蹤期間內的數值積分與求根，可透過替換基線風險函數延伸至其他分布。

| 結果標籤 | 配適模型 | 估計方式 |
|---|---|---|
| exp_PH | Weibull PH 固定形狀為 1 | JM 最大概似 |
| weibull_PH | Weibull PH，估計形狀 | JM 最大概似 |
| spline_PH_lognormal_knots | lognormal 分位數選節點的樣條 PH | JMbayes 貝氏估計 |
| spline_PH_loglogistic_knots | loglogistic 分位數選節點的樣條 PH | JMbayes 貝氏估計 |

後兩項的基線仍為樣條 PH；lognormal、loglogistic 決定節點來源，不是最終基線的參數分布限制。模型比較同時涉及基線結構與估計方法。

## 01｜套件與基線風險函數

MASS 生成相關的隨機效果；nlme 配適縱向混合模型；survival 提供存活資料格式、Cox、參數存活模型與 KM；JM 與 JMbayes 配適聯合模型。

基線風險使用 Weibull 表達式：

$$h_0(t)=\frac{\texttt{shape}}{\texttt{scale}}\left(\frac{t}{\texttt{scale}}\right)^{\texttt{shape}-1}.$$

shape=1 時退化為 exp 基線，scale=100 對應基線 rate 0.01。

```r
library(MASS)
library(nlme)
library(survival)
library(JM)
library(JMbayes)

haz0_weibull <- function(t, shape, scale, eps = 1e-12) {
  tt <- pmax(t, eps)
  (shape / scale) * (tt / scale)^(shape - 1)
}
```

## 02｜有限追蹤期間的事件生成

對病人 i，瞬間風險、累積風險與存活函數為

$$h_i(u)=h_0(u)\exp\{\eta_i(u)\},\qquad
A_i(t)=\int_0^t h_i(u)\,du,\qquad S_i(t)=\exp[-A_i(t)].$$

抽取一次 $E_i\sim\mathrm{Exp}(1)$，等價於 $-\log U_i$，其中 $U_i\sim U(0,1)$。在設限時間 $C_i$ 以前，若累積風險跨過 E，則求解 $A_i(T_i)=E_i$；若 $A_i(C_i)<E_i$，只記錄 $X_i=C_i,\delta_i=0$。

### 02a｜函數介面與輸入

censor_time 是該病人的有限設限時間；t_upper_init 是第一次嘗試的搜尋右端點。函數回傳一個包含 time、event、Ttrue 的清單，不再只回傳事件時間數值。Ttrue 僅對觀察到事件者填入；設限者為 NA，表示未求出其完整事件時間。

```r
cox_exp <- function(id, shape, scale, eta_fun,
                    t_upper_init = 10,
                    t_upper_max = 1e6,
                    eps = 1e-12,
                    rel.tol = 1e-8,
                    subdivisions = 2000L,
                    A_inf = NULL,
                    haz0_fun = NULL) {
  if (is.null(haz0_fun)) {
    haz0_fun <- function(u) haz0_weibull(u, shape, scale)
  }
  if (!is.function(haz0_fun) || !is.function(eta_fun)) {
    stop("haz0_fun 與 eta_fun 必須是函數。")
  }
  if (length(t_upper_init) != 1L || !is.finite(t_upper_init) || t_upper_init <= 0)
    stop("t_upper_init 必須是單一、有限且大於 0 的數值。")
  if (length(censor_time) != 1L || !is.finite(censor_time) || censor_time <= 0)
    stop("censor_time 必須是單一、有限且大於 0 的設限時間。")
  if (length(rel.tol) != 1L || !is.finite(rel.tol) || rel.tol <= 0 ||
      length(subdivisions) != 1L || !is.finite(subdivisions) ||
      subdivisions < 1 || subdivisions != floor(subdivisions))
    stop("rel.tol 須為正數，subdivisions 須為正整數。")
```

### 02b｜累積風險的數值積分

A(t) 計算 0 到 t 的風險面積。

```r
A <- function(t) {
    if (length(t) != 1L || !is.finite(t) || t < 0)
      stop("A(t) 的 t 必須是單一、有限且非負的數值。")
    if (t == 0) return(0)
    val <- integrate(
      function(u) {
        h0 <- haz0_fun(u)
        eta <- eta_fun(u, id)
        if (length(h0) != length(u) || length(eta) != length(u) ||
            any(!is.finite(h0)) || any(h0 < 0) || any(!is.finite(eta)))
          stop("風險函數須對每個積分時間回傳合法值。")
        h <- h0 * exp(eta)
        if (any(!is.finite(h))) stop("個人風險超出數值範圍。")
        h
      },
      lower = 0, upper = t,
      rel.tol = rel.tol, abs.tol = rel.tol,
      subdivisions = subdivisions, stop.on.error = TRUE
    )$value
    if (!is.finite(val) || val < 0)
      stop("病人 ", id, "：累積風險計算失敗。")
    val
  }
```

### 02c｜累積風險門檻與搜尋上限

E 在每位病人的呼叫中只抽一次。f(t)=A(t)-E；fu 儲存目前右端點 tU 的 f 值。fu 小於 0 表示截至 tU 尚未跨過門檻。

搜尋先從 min(t_upper_init,C) 開始。尚未跨過門檻且未達 C 時，上限加倍，最大到 C。例：初始上限 10、C=35，候選上限為 10、20、35；若 C=0.8，第一次就檢查 0.8。加倍是快速擴張範圍的數值策略。

這段沿用上一個 tU 來計算新上限，但每次 A(tU) 都重新積分 0 到 tU，沒有快取前一段面積。

```r
E <- rexp(1)
  if (!is.finite(E) || E <= 0) stop("抽到無效的累積風險門檻。")
  f <- function(t) A(t) - E
  tU <- min(t_upper_init, censor_time)
  fu <- f(tU)
  while (fu < 0) {
    if (tU >= censor_time) {
      return(list(time = censor_time, event = 0L, Ttrue = NA_real_))
    }
    tU <- min(tU * 2, censor_time)
    fu <- f(tU)
  }
```

### 02d｜求根與事件結果

```r
ans <- uniroot(f, interval = c(0, tU), tol = 1e-8, check.conv = TRUE)$root
  if (!is.finite(ans) || ans <= 0 || ans > censor_time)
    stop("病人 ", id, "：有限事件時間計算失敗。")
  list(time = ans, event = 1L, Ttrue = ans)
}
```

## 03｜縱向資料與右設限存活資料

### 03a｜模擬參數與容器函數

simulate_cox_exp 組織整批資料的生成。n 是病人人數；visit_grid 是預定量測時點；beta_fix 是指標的固定截距與斜率；D 是隨機截距、隨機斜率的共變異矩陣；sigma_y 是量測誤差標準差。gamma、alpha 分別是 Z 與潛在指標對風險的係數。

預設 shape=1、scale=100 對應 exp 基線。

```r
simulate_cox_exp <- function(
    n = 5000,
    visit_grid = seq(0, by = 0.2, length.out = 11),
    beta_fix = c(2.0, -0.5),
    D = matrix(c(0.5, 0.01, 0.01, 0.5), 2, 2, byrow = TRUE),
    sigma_y = 0.1, Z_dist = c("binom", "norm"),
    gamma = 1.0, alpha = 2.5,
    shape = 1.0, scale = 100.0,
    admin_tau = 30, cens_rate = 1.6,
    keep_rule = c("<", "<="), haz0_fun = NULL) {
  Z_dist <- match.arg(Z_dist)
  keep_rule <- match.arg(keep_rule)
  if (length(shape) != 1L ||
      !is.finite(shape) ||
      shape != 1) {
    stop("目前 A_inf 的計算公式只適用於 shape = 1。")
  }
  if (length(scale) != 1L ||
      !is.finite(scale) ||
      scale <= 0) {
    stop("scale 必須是有限且大於 0 的數值。")
  }
```

### 03b｜共變數、隨機效果與指標軌跡

Z 預設為成功機率 0.5 的二元共變數。b 的每列是一位病人，b0、b1 分別為隨機截距與斜率。D 對角線 0.5 是變異數，標準差為 $\sqrt{0.5}$；共變異數 0.01 對應相關係數 0.02。

$$m_i(t)=(\beta_0+b_{0i})+(\beta_1+b_{1i})t,\qquad
\eta_i(t)=\gamma Z_i+\alpha m_i(t).$$

m_it 計算潛在指標軌跡；eta_fun 將其連到存活風險。

```r
Z <- if (Z_dist == "binom") rbinom(n, 1, 0.5) else rnorm(n)
  b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
  colnames(b) <- c("b0", "b1")
  m_it <- function(t, i) {
    (beta_fix[1] + b[i, "b0"]) + (beta_fix[2] + b[i, "b1"]) * t
  }
  eta_fun <- function(t, id) gamma * Z[id] + alpha * m_it(t, id)
```

### 03c｜真實事件時間與右設限資料

先根據個人的風險函數生成真實事件時間，再獨立生成設限時間。

在 exp 基線及直線型指標軌跡下，個人風險可寫為：

```math
h_i(t)=\mathrm{base}_i\exp\left(\mathrm{slope}_i t\right).
```

其中，base 是時間零的個人風險率，slope 是線性預測值隨時間變化的斜率。總累積風險為：

```math
A_i(\infty)=
\begin{cases}
\dfrac{\mathrm{base}_i}{-\mathrm{slope}_i}, & \mathrm{slope}_i<0,\\
\infty, & \mathrm{slope}_i\ge 0.
\end{cases}
```

將此上限傳入 `cox_exp`。若抽取的累積風險門檻不低於有限上限，真實事件時間記為 `Inf`；其餘情況使用數值積分與求根取得有限事件時間。

令行政截止時間為 τ，對應程式中的 `admin_tau`。設限時間、觀察時間與事件指標為：

```math
C_i=\min(C_{\mathrm{ind},i},\tau),\qquad
X_i=\min(T_i,C_i),\qquad
\delta_i=I(T_i\le C_i).
```
```r
  # 此處的累積風險上限公式限於 exp 基線
  if (length(shape) != 1L ||
      !is.finite(shape) ||
      shape != 1) {
    stop("目前 A_inf 公式只適用於 shape = 1。")
  }
  if (!is.null(haz0_fun)) {
    stop("目前 A_inf 公式不適用於自訂 haz0_fun。")
  }
  # 先生成每位病人的真實事件時間
  T_true <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    base <- exp(eta_fun(0, i)) / scale
    slope <- alpha * (beta_fix[2] + b[i, "b1"])
    if (!is.finite(base) ||
        base <= 0 ||
        !is.finite(slope)) {
      stop("病人 ", i, "：風險參數計算失敗。")
    }
    if (slope < 0) {
      A_inf <- base / (-slope)
      if (!is.finite(A_inf) || A_inf <= 0) {
        stop("病人 ", i, "：累積風險上限超出數值範圍。")
      }
    } else {
      A_inf <- Inf
    }
    T_true[i] <- cox_exp(
      id = i,
      shape = shape,
      scale = scale,
      eta_fun = eta_fun,
      t_upper_init = admin_tau,
      A_inf = A_inf
    )
  }
  # 再生成獨立設限時間
  Cind <- rexp(n, rate = cens_rate)
  C <- pmin(Cind, admin_tau)

  # 建立觀察時間與事件指標
  time <- pmin(T_true, C)
  event <- as.integer(
    is.finite(T_true) & (T_true <= C)
  )
```

### 03d｜保留追蹤期間的量測

visit_grid 為 0、0.2、…、2，共 11 個預定量測時點。t_i 是第 i 人保留下來的時間，mu 是對應的潛在指標，y_i 加入獨立常態誤差。

long_list 逐人存放量測表；每人觀察時間不同，因此列數可以不同。

```r
long_list <- vector("list", n)
  for (i in seq_len(n)) {
    if (keep_rule == "<") {
      t_i <- visit_grid[visit_grid < time[i]]
    } else {
      t_i <- visit_grid[visit_grid <= time[i]]
    }
    if (!length(t_i)) next
    mu <- m_it(t_i, i)
    y_i <- mu + rnorm(length(t_i), 0, sigma_y)
    long_list[[i]] <- data.frame(id = i, time = t_i, y = y_i, Z = Z[i])
  }
```

### 03e｜合併與回傳

do.call(rbind,long_list) 按列合併各人的量測。long_df 為一人多列，surv_df 為一人一列。long 的 time 是量測時間，surv 的 time 是事件或設限時間。

censoring_rate 由實際 event=0 的比例計算。
```r
long_df <- if (any(!vapply(long_list, is.null, logical(1)))) {
    do.call(rbind, long_list)
  } else {
    data.frame(id = integer(0), time = numeric(0), y = numeric(0), Z = numeric(0))
  }
  surv_df <- data.frame(id = seq_len(n), time = time, event = event,
                        Z = Z, Ttrue = T_true, C = C)
  list(long = long_df, surv = surv_df, censoring_rate = mean(event == 0))
}
```

## 04｜生成一份固定的 exp 資料

執行資料生成。sim 儲存回傳清單，long、surv 取出兩張表。N=5000 對應 nrow(surv)，不是量測表總列數。

```r
set.seed(2025)
sim <- simulate_cox_exp(n = 5000, shape = 1, scale = 100)
long <- sim$long
surv <- sim$surv
print(head(long))
print(head(surv))
print(data.frame(n_subjects = nrow(surv), n_measurements = nrow(long),
                 censoring_rate = sim$censoring_rate))
```

## 05｜初始縱向模型與存活模型

### 05a｜病人對齊

levels_id 建立共同病人順序，factor 指定相同層級，再依病人與量測時間排序。ix=match(long$id,surv$id) 將每筆量測對到該人的存活資料，檢查共變數與追蹤時間一致。

```r
stopifnot(!anyNA(long$id), !anyNA(surv$id), !anyDuplicated(surv$id),
          setequal(long$id, surv$id), !anyDuplicated(long[c("id", "time")]),
          all(is.finite(long$y)), all(is.finite(long$time)),
          all(is.finite(surv$time)), all(surv$time > 0))
levels_id <- sort(unique(surv$id))
long$id <- factor(long$id, levels = levels_id)
surv$id <- factor(surv$id, levels = levels_id)
long <- long[order(long$id, long$time), ]
surv <- surv[order(surv$id), ]
ix <- match(long$id, surv$id)
stopifnot(all(long$Z == surv$Z[ix]), all(long$time <= surv$time[ix]))
```

### 05b｜子模型結構

lmefit 配適固定效果 y~time+Z，以及每人的隨機截距、隨機斜率。公式中的 |id 表示按病人分組；生成端未包含 Z 對指標的直接作用，因此縱向 Z 係數的真值為 0。

```r
lmefit <- nlme::lme(
  y ~ time + Z, random = ~ time | id, data = long,
  na.action = na.fail, control = nlme::lmeControl(opt = "optim")
)
coxfit_weibull <- survival::coxph(
  Surv(time, event) ~ Z, data = surv,
  cluster = id, x = TRUE, model = TRUE, na.action = na.fail
)
coxfit_exp <- coxph(
  Surv(time, event) ~ Z,data = surv,
  cluster = id,x = TRUE,model = TRUE
)
```

## 06｜exp 與 Weibull 聯合模型

weibull-PH-aGH 指定 Weibull PH 及自適應高斯–埃米特積分，用於對未觀察隨機效果的概似積分。

exp 配適透過 scaleWB=1 固定相關 Weibull 參數為 1。Weibull 配適則讓形狀由資料估計。jointfit_exp 與 jointfit_weibull 分別保存兩個結果。

```r
jm_control <- list()
jointfit_exp <- JM::jointModel(
  lmeObject = lmefit, survObject = coxfit_exp,
  timeVar = "time", parameterization = "value",
  method = "weibull-PH-aGH", scaleWB = 1, control = jm_control
)
jointfit_weibull <- JM::jointModel(
  lmeObject = lmefit, survObject = coxfit_weibull,
  timeVar = "time", parameterization = "value",
  method = "weibull-PH-aGH", control = jm_control
)
for (fit in list(jointfit_exp, jointfit_weibull)) {
  if (is.null(fit$convergence) || fit$convergence != 0)
    stop("JM 模型未回報成功收斂。")
}
if (!isTRUE(all.equal(as.numeric(jointfit_exp$coefficients$sigma.t), 1)))
  stop("exp 配適沒有把 Weibull 形狀固定為 1。")
print(summary(jointfit_exp))
print(summary(jointfit_weibull))
```

## 07｜lognormal 節點的樣條 PH

### 07a｜估計分位數與節點

ln0 是只有截距的 lognormal 存活模型，用於挑選節點。mu_hat_ln、sd_hat_ln 是 log 時間的位置與尺度。qlnorm 換算成時間分位數。

knots_ln 僅保留觀察時間範圍內的有限節點，排序並去除重複值。

```r
ln0 <- survival::survreg(Surv(time, event) ~ 1, data = surv, dist = "lognormal")
mu_hat_ln <- as.numeric(coef(ln0))
sd_hat_ln <- as.numeric(ln0$scale)
probs <- seq(0.1, 0.9, length.out = 4)
knots_ln <- qlnorm(probs, meanlog = mu_hat_ln, sdlog = sd_hat_ln)
knots_ln <- sort(unique(knots_ln[is.finite(knots_ln) &
                                  knots_ln > min(surv$time) & knots_ln < max(surv$time)]))
if (!length(knots_ln)) stop("沒有可用的 lognormal 內部節點。")
```

### 07b｜配適樣條聯合模型

baseHaz="regression-splines" 決定最終基線為樣條，knots 指定內部節點，ordSpline=4 表示三次樣條

n.iter 為總迭代數，n.burnin 為前段捨棄次數，n.thin 為抽稀間距。
```r
set.seed(2025)
jm_spline_lnknots <- JMbayes::jointModelBayes(
  lmeObject = lmefit, survObject = coxfit_weibull,
  timeVar = "time", param = "td-value", baseHaz = "regression-splines",
  control = list(knots = knots_ln, ordSpline = 4, ObsTimes.knots = FALSE,
                 n.iter = 12000, n.burnin = 2000, n.thin = 5, seed = 2025)
)
print(summary(jm_spline_lnknots))
```

<a id="loglogistic"></a>
## 08｜loglogistic 節點的樣條 PH

### 08a｜loglogistic 分位數

log 時間採 logistic 的位置–尺度表達式，時間分位數為

$$Q_T(p)=\exp\{\hat\mu+\hat\sigma\,Q_{\mathrm{logistic}}(p)\}.$$

```r
ll0 <- survival::survreg(Surv(time, event) ~ 1, data = surv, dist = "loglogistic")
probs <- seq(0.1, 0.9, length.out = 4)
mu_hat_ll <- as.numeric(coef(ll0))
sd_hat_ll <- as.numeric(ll0$scale)
knots_ll <- exp(mu_hat_ll + sd_hat_ll * qlogis(probs))
knots_ll <- sort(unique(knots_ll[is.finite(knots_ll) &
                                  knots_ll > min(surv$time) & knots_ll < max(surv$time)]))
if (!length(knots_ll)) stop("沒有可用的 loglogistic 內部節點。")
```

### 08b｜樣條模型配適

此分支與 lognormal 分支使用相同種類的聯合模型，差別在節點位置。jm_spline_llknots 表示 loglogistic 選節點的樣條 PH，而非嚴格參數式 Loglogistic-PH。

```r
set.seed(2025)
jm_spline_llknots <- JMbayes::jointModelBayes(
  lmeObject = lmefit, survObject = coxfit_weibull,
  timeVar = "time", param = "td-value", baseHaz = "regression-splines",
  control = list(knots = knots_ll, ordSpline = 4, ObsTimes.knots = FALSE,
                 n.iter = 12000, n.burnin = 2000, n.thin = 5, seed = 2025)
)
print(summary(jm_spline_llknots))
```

## 09｜動態條件存活預測

預測目標為

$$\hat p_i(s,t)=\widehat P(T_i>t\mid T_i>s,\mathcal Y_i(s),Z_i),\qquad t=s+\Delta.$$

s=Tstart 為起點，t=Thoriz 為終點，$\mathcal Y_i(s)$ 為起點以前的指標歷史。模型使用同一份配適結果，隨起點更新個人的條件資訊。

已配適的 joint 模型 + 指定病人在預測起點以前的量測 + 預測起點、預測終點 -> survfitJM() -> 每位病人的條件存活機率 S_hat

### 09a｜依時間與欄名取預測機率

m 是某人的預測摘要表，row 找出對應終點的列號。機率須有限且介於 0 與 1。

```r
extract_survprob <- function(m, Thoriz, id) {
  if (is.data.frame(m)) m <- as.matrix(m)
  if (!is.matrix(m) || !all(c("times", "predSurv") %in% colnames(m)))
    stop("病人 ", id, "：預測摘要缺少 times/predSurv。")
  row <- which(abs(m[, "times"] - Thoriz) < 1e-8)
  if (length(row) != 1L) stop("病人 ", id, "：預測終點不唯一。")
  p <- as.numeric(m[row, "predSurv"])
  if (!is.finite(p) || p < 0 || p > 1) stop("病人 ", id, "：機率不合法。")
  p
}

# JM 用逐人呼叫；JMbayes 保留批次預測形式。
# 固定 simulate=FALSE：兩者都是 plug-in，但參數估計來源不同。
```

### 09b｜準備起點以前的資料

ids 是本次評估的病人編號。di 只取 time≤Tstart 的量測，再依 id 與時間排序。S_hat 預先配置每人的機率並以 id 命名。即使最後量測早於起點，已知存活時間仍使用 Tstart。

```r
get_S_hat <- function(jointfit, long_df, ids, Tstart, Thoriz) {
  stopifnot(length(ids) > 0L, !anyNA(ids), !anyDuplicated(ids),
            Tstart > 0, Thoriz > Tstart)
  di <- long_df[long_df$id %in% ids & long_df$time <= Tstart, , drop = FALSE]
  di <- di[order(di$id, di$time), , drop = FALSE]
  if (!all(ids %in% di$id)) stop("部分病人沒有起點以前的量測。")
  stopifnot(all(is.finite(di$time)), all(is.finite(di$y)), all(is.finite(di$Z)))
  S_hat <- setNames(rep(NA_real_, length(ids)), as.character(ids))
```

### 09c｜JMbayes 的批次預測

inherits 判斷模型類別。JMbayes 分支一次傳入多位病人的量測歷史；sf 為預測物件，sm 為每人摘要清單。依病人名稱取出摘要，避免回傳順序不同而配錯人。

```r
if (inherits(jointfit, "JMbayes")) {
    sf <- JMbayes::survfitJM(jointfit, newdata = di, idVar = "id",
                            last.time = Tstart, survTimes = Thoriz, simulate = FALSE)
    sm <- sf$summaries
    if (is.null(names(sm)) || anyDuplicated(names(sm)) ||
        !setequal(names(sm), as.character(ids)))
      stop("JMbayes 回傳的病人編號與輸入不一致。")
    for (k in seq_along(ids)) {
      i <- as.character(ids[k])
      S_hat[k] <- extract_survprob(sm[[i]], Thoriz, i)
    }
```

### 09d｜JM 的逐人預測與回傳

JM 分支逐人呼叫 survfitJM。此設定在估計參數下使用個人隨機效果的經驗貝氏後驗眾數進行代入式預測。

函數最後回傳 S_hat。

```r
} else if (inherits(jointfit, "jointModel")) {
    for (k in seq_along(ids)) {
      i <- ids[k]
      sf <- JM::survfitJM(jointfit, newdata = di[di$id == i, , drop = FALSE],
                          idVar = "id", last.time = Tstart,
                          survTimes = Thoriz, simulate = FALSE)
      S_hat[k] <- extract_survprob(sf$summaries[[1]], Thoriz, i)
    }
  } else {
    stop("目前預測介面僅支援 JM 與 JMbayes，不含 JMbayes2/JSM。")
  }
  S_hat
}
```

## 10｜單一起點的 IPCW Brier Score

### 10a｜風險集與機率對齊

只納入 $X_i>s$ 且起點以前有量測資料的人，形成 $\mathcal R_s$。R 是其存活資料表，n_risk 是人數 $n_s$。

S_hat 依 R 的 id 順序重新對齊。預測時不使用起點以後的指標，觀察結局則用於後續評分。

```r
# 一個起點、一個終點的評分：條件反向 KM +  IPCW 貢獻。
brier_ipcw <- function(jointfit, long_df, surv_df, Tstart, Thoriz, min_G = 0.05) {
  stopifnot(length(Tstart) == 1L, is.finite(Tstart), Tstart > 0,
            length(Thoriz) == 1L, is.finite(Thoriz), Thoriz > Tstart,
            length(min_G) == 1L, is.finite(min_G), min_G > 0, min_G <= 1,
            !anyNA(surv_df$id), !anyDuplicated(surv_df$id),
            all(is.finite(surv_df$time)), all(surv_df$time > 0),
            all(surv_df$event %in% c(0, 1)))
  ids <- surv_df$id[surv_df$time > Tstart]
  available <- unique(long_df$id[long_df$time <= Tstart])
  n_no_history <- sum(!ids %in% available)
  ids <- ids[ids %in% available]
  if (length(ids) < 2L) stop("當次評估病人不足 2 人。")
  R <- surv_df[match(ids, surv_df$id), , drop = FALSE]
  S_hat <- get_S_hat(jointfit, long_df, ids, Tstart, Thoriz)
  S_hat <- unname(S_hat[as.character(R$id)])
  stopifnot(length(S_hat) == nrow(R), all(is.finite(S_hat)),
            all(S_hat >= 0 & S_hat <= 1))
```

### 10b｜三種結局與反向 KM

令 $X_i=\min(T_i,C_i)$、$\delta_i=I(T_i\le C_i)$。在起點風險集中，定義 $G_s(u)=P(C_i>u\mid C_i>s)$。

| 條件 | 程式名稱 | 已知結局 $I(T_i>t)$ | IPCW 權重 |
|---|---|---:|---:|
| $s<X_i\le t,\delta_i=1$ | dead | 0 | $1/\hat G_s(X_i-)$ |
| $X_i>t$ | alive | 1 | $1/\hat G_s(t)$ |
| $s<X_i\le t,\delta_i=0$ | unknown | 未知 | 0 |

Surv(time,1-event) 以設限作為 KM 的事件、真正事件作為 KM 的設限，估計設限存活機率。

G_step 依 KM 跳點查值：left=TRUE 使用嚴格小於，得到事件前一瞬間的左極限；一般查終點用小於等於。在相鄰跳點之間保持上一個 KM 值。

```r
dead <- R$event == 1 & R$time <= Thoriz
  alive <- R$time > Thoriz
  unknown <- !dead & !alive
  km <- survival::survfit(survival::Surv(time, 1 - event) ~ 1, data = R)
  G_step <- function(u, left = FALSE) {
    vapply(u, function(x) {
      idx <- if (left) which(km$time < x) else which(km$time <= x)
      if (length(idx)) tail(km$surv[idx], 1) else 1
    }, numeric(1))
  }
```

### 10c｜加權平方誤差

$$\widehat{BS}(s,t)=\frac1{n_s}\sum_{i\in\mathcal R_s}
\left[
\frac{I(X_i\le t,\delta_i=1)\hat p_i(s,t)^2}{\hat G_s(X_i-)}
+\frac{I(X_i>t)\{1-\hat p_i(s,t)\}^2}{\hat G_s(t)}
\right].$$

contrib 儲存每人的貢獻，unknown 留為 0，但仍留在分母 n_risk。分母不是已知結局的人數，也不是權重總和。min_G=0.05 是本實作的支持度停止門檻，不是將較小的 G 裁切到 0.05。

```r
G_ti <- G_step(R$time[dead], left = TRUE)
  G_t <- G_step(Thoriz)
  denom <- c(G_ti, if (any(alive)) G_t else numeric(0))
  if (length(denom) && any(!is.finite(denom) | denom < min_G))
    stop("設限權重缺乏足夠支持；不可裁切機率或刪人後繼續。")
  contrib <- numeric(nrow(R))
  contrib[dead] <- S_hat[dead]^2 / G_ti
  contrib[alive] <- (1 - S_hat[alive])^2 / G_t
  data.frame(Tstart = Tstart, Thoriz = Thoriz, Dt = Thoriz - Tstart,
             Brier = mean(contrib), n_risk = nrow(R),
             n_event = sum(dead), n_alive = sum(alive), n_unknown = sum(unknown),
             n_no_history = n_no_history,
             min_G = if (length(denom)) min(denom) else NA_real_)
}
```

## 11｜移動起點、固定跨度

brier_grid 對每個 Tstarts 執行同一個評分步驟。j 為起點索引，s 是當次起點，t=s+Dt 是終點；out[[j]] 保存一列結果，最後以 rbind 合併。

每個起點重新選人、讀取歷史、計算條件預測與 KM 權重。不重新配適 joint 模型，不沿用前次預測機率，也不累加前次 Brier。tryCatch 補上出錯起點與終點後停止，避免默默略過失敗的時間。

```r
# 移動 Tstart，固定 Dt；不重新配適模型。
brier_grid <- function(jointfit, long_df, surv_df, Tstarts, Dt, min_G = 0.05) {
  stopifnot(length(Tstarts) > 0L, all(is.finite(Tstarts)), all(Tstarts > 0),
            all(diff(Tstarts) > 0), length(Dt) == 1L, is.finite(Dt), Dt > 0)
  out <- vector("list", length(Tstarts))
  for (j in seq_along(Tstarts)) {
    s <- Tstarts[j]
    t <- s + Dt
    message("起點 ", j, "/", length(Tstarts), "：", s, " -> ", t)
    out[[j]] <- tryCatch(
      brier_ipcw(jointfit, long_df, surv_df, s, t, min_G),
      error = function(e) stop("Tstart=", s, ", Thoriz=", t, "：",
                               conditionMessage(e), call. = FALSE)
    )
  }
  do.call(rbind, out)
}
```

## 12｜四個模型的完整時間網格與輸出

起點由 0.0001 開始，每次增加 0.2，共 12 個，最後起點為 2.2001，終點為 2.4001。近零起點不完全等於 0。

```r
Tstarts <- seq(0.0001, 2.4, by = 0.2)
Dt <- 0.2
fits <- list(
  exp_PH = jointfit_exp,
  weibull_PH = jointfit_weibull,
  spline_PH_lognormal_knots = jm_spline_lnknots,
  spline_PH_loglogistic_knots = jm_spline_llknots
)
tables <- lapply(fits, function(fit) brier_grid(fit, long, surv, Tstarts, Dt))
tab_exp <- tables$exp_PH
tab_weibull <- tables$weibull_PH
tab_spline_lnknots <- tables$spline_PH_lognormal_knots
tab_spline_llknots <- tables$spline_PH_loglogistic_knots
comparison <- do.call(rbind, lapply(names(tables), function(name) {
  cbind(model = name, tables[[name]])
}))
```

## 13｜半參數聯合模型的目的與 start–stop 資料

### 13.1｜分析目的

本節使用同一份 exp 基線生成資料，增加半參數 Cox 聯合模型作為比較。

模型仍保留潛在指標與事件風險的關聯：

$$
h_i(t)=h_0(t)\exp\{\gamma Z_i+\alpha m_i(t)\}.
$$

與 exp、Weibull 配適的主要差別，在於半參數模型不指定基線風險必須符合某個參數分布。迴歸係數與指標關聯係數仍需估計，因此稱為「半參數」。

本節使用 JSM 套件的 jmodelTM。設定 model=1，以當下潛在指標值連結事件風險；rho=0 選擇 Cox PH。

分析流程為：

1. 將量測與存活資料整理為 start–stop 格式。
2. 建立縱向及 Cox 初始模型。
3. 配適 JSM 半參數聯合模型。
4. 使用起點以前的歷史資料進行動態預測。
5. 計算固定預測跨度的 IPCW Brier score。

### 13.2｜start–stop 的資料意義

start–stop 將一位病人的追蹤期間分成多列。

例如病人在 0、0.2、0.4 有量測，於 0.55 發生事件：

| start | stop | y 的量測時間 | event |
|---:|---:|---:|---:|
| 0 | 0.2 | 0 | 0 |
| 0.2 | 0.4 | 0.2 | 0 |
| 0.4 | 0.55 | 0.4 | 1 |

start 是區段起點，stop 是區段終點，event 表示是否在該段終點發生事件。若最後在 0.55 設限，最後一列的 event 也為 0。

這仍是右設限資料。Surv(start, stop, event) 表示追蹤區段，不表示只知道事件發生在兩個端點之間。

JSM 使用相同資料表配適縱向與存活子模型，因此需要這種資料排列。每列保存一筆量測，但最終潛在指標仍由混合模型的截距與斜率描述，不是把量測值當作整段不變的真實指標。

### 13.3｜連結量測與存活資料

long 為一人多列的量測表，surv 為一人一列的存活表。

建立副本後，將每筆量測對應到該病人的追蹤終點。Stime 保存追蹤終點，status 保存病人層級的事件指標。

此處沿用量測時間嚴格早於追蹤終點的設定。

```r
library(nlme)
library(survival)
library(JSM)

long_jsm <- long[, c("id", "time", "y", "Z")]
surv_jsm <- surv[, c("id", "time", "event", "Z")]

long_jsm$id <- as.character(long_jsm$id)
surv_jsm$id <- as.character(surv_jsm$id)

stopifnot(
  !anyNA(long_jsm$id),
  !anyNA(surv_jsm$id),
  !anyDuplicated(surv_jsm$id),
  setequal(long_jsm$id, surv_jsm$id),
  !anyDuplicated(long_jsm[c("id", "time")]),
  all(is.finite(long_jsm$time)),
  all(is.finite(long_jsm$y)),
  all(is.finite(long_jsm$Z)),
  all(is.finite(surv_jsm$time)),
  all(is.finite(surv_jsm$Z)),
  all(surv_jsm$time > 0),
  all(surv_jsm$event %in% c(0, 1))
)

ix <- match(long_jsm$id, surv_jsm$id)

stopifnot(all(long_jsm$Z == surv_jsm$Z[ix]))

long_jsm$Stime <- surv_jsm$time[ix]
long_jsm$status <- surv_jsm$event[ix]

if (any(long_jsm$time < 0 |
        long_jsm$time >= long_jsm$Stime)) {
  stop("量測時間必須非負，且嚴格早於追蹤終點。")
}
```

### 13.4｜逐人建立追蹤區段

di 是單一病人的量測資料，n_i 是其量測筆數。

每列的 start 為量測時間，stop 為下一次量測時間；最後一列延伸到 Stime。最後一段不新增量測，也不複製額外的 y。

event_i 先全部設為 0，再將最後一列設為該病人的事件狀態。

```r
ids <- sort(unique(surv_jsm$id))
interval_list <- vector("list", length(ids))

for (k in seq_along(ids)) {

  i <- ids[k]

  di <- long_jsm[long_jsm$id == i, , drop = FALSE]
  di <- di[order(di$time), , drop = FALSE]

  if (di$time[1] != 0) {
    stop("病人 ", i, "：第一筆量測須位於時間 0。")
  }

  n_i <- nrow(di)
  # 建立 區間終點
  stop_i <- c(di$time[-1], di$Stime[1])

  event_i <- integer(n_i)
  event_i[n_i] <- di$status[1]

  interval_list[[k]] <- data.frame(
    id = i,
    start = di$time,
    stop = stop_i,
    y = di$y,
    Z = di$Z,
    Stime = di$Stime,
    event = event_i
  )
}
```

### 13.5｜合併區段資料

data_long 是 JSM 配適使用的資料表。id 設為具有固定順序的 factor，並依病人與區段起點排序。

整理後，每筆原始量測仍對應一列，事件總筆數也應保持一致。

```r
data_long <- do.call(rbind, interval_list)
rownames(data_long) <- NULL

data_long$id <- factor(data_long$id, levels = ids)

data_long <- data_long[
  order(data_long$id, data_long$start),
  ,
  drop = FALSE
]

stopifnot(
  all(data_long$start < data_long$stop),
  nrow(data_long) == nrow(long_jsm),
  sum(data_long$event) == sum(surv_jsm$event)
)

head(data_long, 10)
```

## 14｜配適 JSM 半參數聯合模型

### 14.1｜建立初始子模型

縱向初始模型包含固定截距、時間及 Z，以及個人的隨機截距與斜率。

依 JSM 文件，輸入的 lme 使用 pdDiag 結構。

Cox 初始模型使用 start–stop 格式。x=TRUE 保存設計矩陣，y=TRUE 保存實際使用的存活反應資料。

```r
lme_model <- nlme::lme(
  fixed = y ~ start + Z,
  random = list(id = nlme::pdDiag(~ start)),
  data = data_long,
  na.action = na.fail,
  control = nlme::lmeControl(opt = "optim")
)

cox_model <- survival::coxph(
  survival::Surv(start, stop, event) ~ Z,
  data = data_long,
  x = TRUE,
  y = TRUE,
  model = TRUE,
  na.action = na.fail
)
```

### 14.2｜聯合模型與積分節點

jmodelTM 共同估計縱向與存活部分。

nknot 控制每個隨機效果維度的高斯–埃米特積分節點數。

此處使用較少節點以降低計算負擔，因此將結果命名為 jm_model_cox_k3。

```r
gc()

jm_model_cox_k3 <- jmodelTM(
  fitLME = lme_model,
  fitCOX = cox_model,
  data = data_long,
  model = 1,
  rho = 0,
  timeVarY = "start",
  control = list(
    nknot = 3,
    max.iter = 250,
    tol.P = 1e-3,
    tol.L = 1e-6
  )
)

print(jm_model_cox_k3$convergence)

if (!identical(jm_model_cox_k3$convergence, "success")) {
  stop("JSM 未回報成功收斂，暫不進行預測。")
}

summary(jm_model_cox_k3)
```

## 15｜整理聯合模型的預測參數

### 15.1｜參數與預測範圍

後續預測使用最終 joint model 的估計，而非初始 lme_model 的估計。

| 代稱 | 意義 |
|---|---|
| beta | 縱向固定截距、時間與 Z 係數 |
| phi | Z 對事件風險的係數 |
| alpha | 潛在指標對事件風險的關聯係數 |
| sigma | 指標量測誤差標準差 |
| D | 隨機截距與斜率的共變異矩陣 |
| u | 基線累積風險跳躍的事件時間 |
| dH0 | 各事件時間的基線累積風險增量 |

本節預測程式限定目前的模型結構：

- 縱向：y ~ start + Z。
- 隨機效果：截距與 start 斜率。
- 存活共變數：Z。
- 關聯方式：model=1。
- 存活模型：rho=0。

JSM 的基線增量直接加總，不再乘上時間間距。

### 15.2｜取出並檢查參數

```r
prepare_jsm_prediction <- function(jsm_fit) {

  if (!identical(jsm_fit$convergence, "success")) {
    stop("JSM 尚未成功收斂。")
  }

  if (length(jsm_fit$rho) != 1L || jsm_fit$rho != 0) {
    stop("此程式只適用於 rho = 0。")
  }

  cc <- jsm_fit$coefficients
  beta_names <- c("(Intercept)", "start", "Z")

  if (!setequal(names(cc$beta), beta_names) ||
      length(cc$phi) != 1L ||
      !identical(names(cc$phi), "Z")) {
    stop("模型固定效果結構與目前預測程式不一致。")
  }

  D <- cc[["Bsigma"]]
  if (is.null(D)) D <- cc[["BSigma"]]

  re_names <- c("(Intercept)", "start")

  if (!is.matrix(D) ||
      !all(dim(D) == c(2, 2)) ||
      !all(re_names %in% rownames(D)) ||
      !all(re_names %in% colnames(D))) {
    stop("無法確認隨機截距與斜率的共變異矩陣。")
  }

  D <- D[re_names, re_names, drop = FALSE]

  lamb <- as.data.frame(cc$lamb)

  if (!all(c("time", "bashaz") %in% names(lamb))) {
    stop("基線估計欄位與預期不同。")
  }

  lamb <- lamb[order(lamb$time), , drop = FALSE]

  p <- list(
    beta = cc$beta[beta_names],
    phi = as.numeric(cc$phi),
    alpha = as.numeric(cc$alpha),
    sigma = as.numeric(cc$Ysigma),
    D = D,
    u = lamb$time,
    dH0 = lamb$bashaz
  )

  stopifnot(
    all(is.finite(p$beta)),
    is.finite(p$phi),
    length(p$alpha) == 1L,
    is.finite(p$alpha),
    length(p$sigma) == 1L,
    is.finite(p$sigma),
    p$sigma > 0,
    all(is.finite(D)),
    isSymmetric(D),
    min(eigen(D, symmetric = TRUE)$values) > 0,
    length(p$u) > 0L,
    all(is.finite(p$u)),
    all(p$u > 0),
    all(diff(p$u) > 0),
    all(is.finite(p$dH0)),
    all(p$dH0 >= 0)
  )

  p$Dinv <- solve(D)
  p
}

jsm_par <- prepare_jsm_prediction(jm_model_cox_k3)
```

## 16｜個人累積風險與動態預測

### 16.1｜未來窗口的累積風險

給定個人的隨機效果 b，潛在指標為：

```math
m_i(u)=
\widehat{\beta}_0
+\widehat{\beta}_1 u
+\widehat{\beta}_Z Z_i
+b_{0i}
+b_{1i}u.
```

在起點 s 到終點 t 的累積風險增量為：

```math
H_i(s,t\mid b_i)=\sum_{s \lt u_j \le t}\Delta\widehat{H}_0(u_j)\exp\left(\widehat{\phi}Z_i+\widehat{\alpha}m_i(u_j)\right).
```

其中，u 的下標 j 表示基線累積風險的跳點；每個跳點的基線累積風險增量，乘上個人的風險倍率，再於指定窗口內加總。

`jsm_cumhaz` 的 `b` 是兩個個人隨機效果，`z` 是該人的 Z，`p` 是前段整理的模型參數。`lower` 與 `upper` 分別是窗口起點與終點。

```r
jsm_cumhaz <- function(b, z, lower, upper, p) {

  use <- p$u > lower &
         p$u <= upper &
         p$dH0 > 0

  if (!any(use)) return(0)

  u <- p$u[use]

  m <- p$beta["(Intercept)"] +
       p$beta["start"] * u +
       p$beta["Z"] * z +
       b[1] + b[2] * u

  log_terms <- log(p$dH0[use]) +
               p$phi * z +
               p$alpha * m

  a <- max(log_terms)

  if (!is.finite(a)) {
    stop("累積風險計算出現非有限值。")
  }

  log_total <- a + log(sum(exp(log_terms - a)))

  if (log_total > log(.Machine$double.xmax)) {
    return(Inf)
  }

  exp(log_total)
}
```
### 16.2｜由起點以前的資料估計個人效果

預測目標為：

$$
\hat p_i(s,t)=
\widehat P(T_i>t\mid T_i>s,Y_i(\le s),Z_i).
$$

本節採用經驗貝氏後驗眾數代入式預測。個人效果依下列分布估計：

$$
p(b_i\mid Y_i(\le s),T_i>s)
\propto
p(Y_i(\le s)\mid b_i)\,p(b_i)\,S_i(s\mid b_i).
$$

因此不能直接使用配適物件中以完整個人資料估計的 est.bi。

以下函數先截取起點以前的量測。Q 為隨機效果設計矩陣，b_long、V 是只考慮這些量測時的常態後驗中心與共變異矩陣。

```r
predict_S_jsm <- function(di, z, Tstart, Thoriz, p) {

  di <- di[di$time <= Tstart, , drop = FALSE]
  di <- di[order(di$time), , drop = FALSE]

  if (nrow(di) == 0L) {
    stop("起點以前沒有量測資料。")
  }

  stopifnot(
    Thoriz > Tstart,
    all(is.finite(di$time)),
    all(is.finite(di$y)),
    all(di$Z == z)
  )

  Q <- cbind(1, di$time)

  fixed_mean <- p$beta["(Intercept)"] +
                p$beta["start"] * di$time +
                p$beta["Z"] * z

  residual <- di$y - fixed_mean

  V <- solve(
    p$Dinv + crossprod(Q) / p$sigma^2
  )

  b_long <- as.vector(
    V %*% crossprod(Q, residual) / p$sigma^2
  )

  L <- t(chol(V))
```

### 16.3｜加入存活到起點的資訊

以標準化座標 a 表示個人效果：

$$b=b_{\mathrm{long}}+La.$$

objective 是負對數後驗密度中與 a 有關的部分。第一項來自量測後驗，第二項為存活到起點的累積風險。

```r
  objective <- function(a) {

    b <- b_long + as.vector(L %*% a)

    H_start <- jsm_cumhaz(
      b = b,
      z = z,
      lower = 0,
      upper = Tstart,
      p = p
    )

    sum(a^2) / 2 + H_start
  }

  if (!is.finite(objective(c(0, 0)))) {
    stop("個人效果的最佳化起始值無效。")
  }

  opt <- optim(
    par = c(0, 0),
    fn = objective,
    method = "BFGS",
    control = list(
      maxit = 500,
      reltol = 1e-10
    )
  )

  if (opt$convergence != 0L ||
      !is.finite(opt$value) ||
      any(!is.finite(opt$par))) {
    stop("個人隨機效果估計未收斂。")
  }

  b_hat <- b_long + as.vector(L %*% opt$par)
```

### 16.4｜取得條件存活機率

將 b_hat 代入累積風險，得到：

$$
\hat S_i(t\mid s)=\exp\{-H_i(s,t\mid\hat b_i(s))\}.
$$

S_hat 是一位病人在該窗口的預測機率。

```r
  H_future <- jsm_cumhaz(
    b = b_hat,
    z = z,
    lower = Tstart,
    upper = Thoriz,
    p = p
  )

  S_hat <- exp(-H_future)

  if (!is.finite(S_hat) || S_hat < 0 || S_hat > 1) {
    stop("預測機率不合法。")
  }

  S_hat
}
```

## 17｜JSM 的 IPCW Brier Score

### 17.1｜風險集與逐人預測

起點仍被追蹤，且起點以前有量測的人。R 是當次評估資料，history 是依病人分組的量測清單。

n_no_history 記錄缺少量測歷史而未納入的人。

```r
brier_jsm_ipcw <- function(
    p, history, surv_df,
    Tstart, Thoriz, min_G = 0.05) {

  stopifnot(
    length(Tstart) == 1L,
    is.finite(Tstart),
    Tstart > 0,
    length(Thoriz) == 1L,
    is.finite(Thoriz),
    Thoriz > Tstart,
    length(min_G) == 1L,
    is.finite(min_G),
    min_G > 0,
    min_G <= 1
  )

  R <- surv_df[
    surv_df$time > Tstart,
    ,
    drop = FALSE
  ]

  has_history <- vapply(
    as.character(R$id),
    function(i) {
      di <- history[[i]]
      !is.null(di) && any(di$time <= Tstart)
    },
    logical(1)
  )

  n_no_history <- sum(!has_history)
  R <- R[has_history, , drop = FALSE]

  if (nrow(R) < 2L) {
    stop("起點符合條件的病人不足 2 人。")
  }

  S_hat <- vapply(seq_len(nrow(R)), function(k) {

    i <- as.character(R$id[k])

    tryCatch(
      predict_S_jsm(
        di = history[[i]],
        z = R$Z[k],
        Tstart = Tstart,
        Thoriz = Thoriz,
        p = p
      ),
      error = function(e) {
        stop(
          "病人 ", i, "：", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }, numeric(1))
```

### 17.2｜三種情況與設限存活函數

令 X 為觀察時間，δ 為事件指標，條件設限存活函數為：

$$G_s(u)=P(C>u\mid C>s).$$

| 情況 | 名稱 | 評分貢獻 |
|---|---|---|
| 窗口內發生事件 | dead | $\hat S_i^2/\hat G_s(X_i-)$ |
| 追蹤超過終點 | alive | $(1-\hat S_i)^2/\hat G_s(t)$ |
| 窗口內設限 | unknown | 0，仍保留在平均分母 |

Surv(time,1-event) 將設限當作 KM 事件，用反向 KM 估計設限存活機率。

本分析在當次風險集上估計條件 G_s，因此不再額外乘 G(s)。這依賴目前獨立設限的設定。

G_step 的 left=TRUE 取得事件前一瞬間的左極限；查終點時則使用該時間的階梯函數值。

```r
  dead <- R$event == 1 & R$time <= Thoriz
  alive <- R$time > Thoriz
  unknown <- !dead & !alive

  km <- survival::survfit(
    survival::Surv(time, 1 - event) ~ 1,
    data = R
  )

  G_step <- function(x, left = FALSE) {

    vapply(x, function(u) {

      index <- if (left) {
        which(km$time < u)
      } else {
        which(km$time <= u)
      }

      if (length(index)) {
        tail(km$surv[index], 1)
      } else {
        1
      }
    }, numeric(1))
  }

  G_event <- G_step(R$time[dead], left = TRUE)
  G_end <- G_step(Thoriz)

  denominators <- c(
    G_event,
    if (any(alive)) G_end else numeric(0)
  )

  if (any(!is.finite(denominators) |
          denominators < min_G)) {
    stop("設限存活機率過小，該窗口缺乏足夠評分支持。")
  }
```

### 17.3｜平均加權平方誤差

contribution 保存每人的加權誤差。未知結局者貢獻為 0，但仍包含在 n_risk 中。

min_G 是程式的支持度停止門檻。

```r
  contribution <- numeric(nrow(R))

  contribution[dead] <- S_hat[dead]^2 / G_event

  if (any(alive)) {
    contribution[alive] <-
      (1 - S_hat[alive])^2 / G_end
  }

  data.frame(
    Tstart = Tstart,
    Thoriz = Thoriz,
    Dt = Thoriz - Tstart,
    Brier = mean(contribution),
    n_risk = nrow(R),
    n_event = sum(dead),
    n_alive = sum(alive),
    n_unknown = sum(unknown),
    n_no_history = n_no_history
  )
}
```

## 18｜執行完整時間網格與模型比較

### 18.1｜整理預測使用的量測歷史

history_jsm 將量測依病人分組，方便逐人讀取。雖然清單中保存完整量測，predict_S_jsm 每次仍只使用 Tstart 以前的資料。

```r
stopifnot(
  !anyNA(surv$id),
  !anyDuplicated(surv$id),
  !anyNA(long$id),
  !anyDuplicated(long[c("id", "time")]),
  all(is.finite(surv$time)),
  all(surv$time > 0),
  all(surv$event %in% c(0, 1)),
  all(is.finite(surv$Z)),
  all(is.finite(long$time)),
  all(long$time >= 0),
  all(is.finite(long$y)),
  all(is.finite(long$Z))
)

history_jsm <- split(
  long,
  as.character(long$id),
  drop = TRUE
)

Tstarts <- seq(0.0001, 2.4, by = 0.2)
Dt <- 0.2

if (max(Tstarts + Dt) > max(surv$time)) {
  stop("預測終點超過目前資料的最長追蹤時間。")
}
```

### 18.2｜移動起點、固定跨度

每次起點向後移動 0.2，終點為起點加上 Dt。每個窗口重新選擇風險集、估計個人效果、預測與評分，不重新配適整個聯合模型。

```r
result_jsm <- vector("list", length(Tstarts))

for (j in seq_along(Tstarts)) {

  s <- Tstarts[j]
  t <- s + Dt

  message(
    "正在計算第 ", j, "/", length(Tstarts),
    " 組：", s, " → ", t
  )

  result_jsm[[j]] <- tryCatch(
    brier_jsm_ipcw(
      p = jsm_par,
      history = history_jsm,
      surv_df = surv,
      Tstart = s,
      Thoriz = t
    ),
    error = function(e) {
      stop(
        "Tstart = ", s,
        "，Thoriz = ", t,
        "：", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

tab_semi_k3 <- do.call(rbind, result_jsm)

print(
  tab_semi_k3,
  digits = 10,
  row.names = FALSE
)
```

### 18.3｜結果解讀與比較範圍

若完整執行目前設定的 12 組起點與終點，`tab_semi_k3` 成功完成後共有 12 列。每列對應一組起點與終點，而非單一病人。

- `Brier`：該窗口的 IPCW 平均平方預測誤差。
- `n_risk`：納入評分的病人人數。
- `n_event`：窗口內發生事件的人數。
- `n_alive`：追蹤超過終點的人數。
- `n_unknown`：窗口內設限、終點結局未知的人數。
- `n_no_history`：因缺少起點以前量測而未納入的人數。

各列應滿足：

```math
n_{\mathrm{event}}
+
n_{\mathrm{alive}}
+
n_{\mathrm{unknown}}
=
n_{\mathrm{risk}}.
```

與 exp、Weibull、樣條 PH 比較時，需使用同一份資料、相同時間網格、相同風險集選取方式與 IPCW 規則。

目前結果具有以下範圍：

1. 使用配適資料本身評分，屬於樣本內表現。
2. 個人效果使用起點以前歷史與存活資訊重新估計，但整體模型參數仍由完整訓練資料估計。
3. 預測為經驗貝氏後驗眾數代入，沒有完整積分個人效果及模型參數的不確定性。
4. JSM 配適採 `nknot = 3`，收斂不等於積分精度已充分確認。
5. 本節為自訂 JSM 預測與評分實作，不能僅以程式成功執行視為數值及統計方法已驗證。
6. 此處尚未計算 IBS；沿移動起點積分與固定起點、沿終點積分是不同設計。

# 結果
## Exp

| Tstart<dbl> | Thoriz<dbl> | Dt<dbl> | Brier<dbl> | n_risk<int> | n_event<int> | n_alive<int> | n_unknown<int> |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.0001 | 0.2001 | 0.2 | 0.132191246 | 4993 | 2010 | 2013 | 970 |
| 0.2001 | 0.4001 | 0.2 | 0.114821440 | 2013 | 326 | 1162 | 525 |
| 0.4001 | 0.6001 | 0.2 | 0.081514825 | 1162 | 118 | 744 | 300 |
| 0.6001 | 0.8001 | 0.2 | 0.059226997 | 744 | 56 | 488 | 200 |
| 0.8001 | 1.0001 | 0.2 | 0.047423178 | 488 | 23 | 338 | 127 |
| 1.0001 | 1.2001 | 0.2 | 0.054885274 | 338 | 21 | 239 | 78 |
| 1.2001 | 1.4001 | 0.2 | 0.015032535 | 239 | 5 | 166 | 68 |
| 1.4001 | 1.6001 | 0.2 | 0.015514769 | 166 | 3 | 125 | 38 |
| 1.6001 | 1.8001 | 0.2 | 0.043579271 | 125 | 5 | 83 | 37 |
| 1.8001 | 2.0001 | 0.2 | 0.012163100 | 83 | 1 | 61 | 21 |
| 2.0001 | 2.2001 | 0.2 | 0.015436711 | 61 | 1 | 47 | 13 |
| 2.2001 | 2.4001 | 0.2 | 0.001752082 | 47 | 0 | 37 | 10 |

## Weibull

| Tstart<dbl> | Thoriz<dbl> | Dt<dbl> | Brier<dbl> | n_risk<int> | n_event<int> | n_alive<int> | n_unknown<int> |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.0001 | 0.2001 | 0.2 | 0.13218727 | 4993 | 2010 | 2013 | 970 |
| 0.2001 | 0.4001 | 0.2 | 0.11482543 | 2013 | 326 | 1162 | 525 |
| 0.4001 | 0.6001 | 0.2 | 0.08151880 | 1162 | 118 | 744 | 300 |
| 0.6001 | 0.8001 | 0.2 | 0.05924968 | 744 | 56 | 488 | 200 |
| 0.8001 | 1.0001 | 0.2 | 0.04740874 | 488 | 23 | 338 | 127 |
| 1.0001 | 1.2001 | 0.2 | 0.05490367 | 338 | 21 | 239 | 78 |
| 1.2001 | 1.4001 | 0.2 | 0.01504612 | 239 | 5 | 166 | 68 |
| 1.4001 | 1.6001 | 0.2 | 0.01551694 | 166 | 3 | 125 | 38 |
| 1.6001 | 1.8001 | 0.2 | 0.04356819 | 125 | 5 | 83 | 37 |
| 1.8001 | 2.0001 | 0.2 | 0.01216580 | 83 | 1 | 61 | 21 |
| 2.0001 | 2.2001 | 0.2 | 0.01543527 | 61 | 1 | 47 | 13 |
| 2.2001 | 2.4001 | 0.2 | 0.00173658 | 47 | 0 | 37 | 10 |

## lognormal


| Tstart | Thoriz | Dt | Brier | n_risk | n_event | n_alive | n_unknown |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0.0001 | 0.2001 | 0.2 | 0.1323299376 | 4993 | 2010 | 2013 | 970 |
| 0.2001 | 0.4001 | 0.2 | 0.1148031737 | 2013 | 326 | 1162 | 525 |
| 0.4001 | 0.6001 | 0.2 | 0.0815246180 | 1162 | 118 | 744 | 300 |
| 0.6001 | 0.8001 | 0.2 | 0.0591012628 | 744 | 56 | 488 | 200 |
| 0.8001 | 1.0001 | 0.2 | 0.0474421610 | 488 | 23 | 338 | 127 |
| 1.0001 | 1.2001 | 0.2 | 0.0549005760 | 338 | 21 | 239 | 78 |
| 1.2001 | 1.4001 | 0.2 | 0.0149896225 | 239 | 5 | 166 | 68 |
| 1.4001 | 1.6001 | 0.2 | 0.0155222703 | 166 | 3 | 125 | 38 |
| 1.6001 | 1.8001 | 0.2 | 0.0437410797 | 125 | 5 | 83 | 37 |
| 1.8001 | 2.0001 | 0.2 | 0.0121035655 | 83 | 1 | 61 | 21 |
| 2.0001 | 2.2001 | 0.2 | 0.0154376753 | 61 | 1 | 47 | 13 |
| 2.2001 | 2.4001 | 0.2 | 0.0030924884 | 47 | 0 | 37 | 10 |


## loglogistic


| Tstart | Thoriz | Dt | Brier | n_risk | n_event | n_alive | n_unknown |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0.0001 | 0.2001 | 0.2 | 0.132273234023 | 4993 | 2010 | 2013 | 970 |
| 0.2001 | 0.4001 | 0.2 | 0.114821232817 | 2013 | 326 | 1162 | 525 |
| 0.4001 | 0.6001 | 0.2 | 0.081562733563 | 1162 | 118 | 744 | 300 |
| 0.6001 | 0.8001 | 0.2 | 0.059365120356 | 744 | 56 | 488 | 200 |
| 0.8001 | 1.0001 | 0.2 | 0.047219009528 | 488 | 23 | 338 | 127 |
| 1.0001 | 1.2001 | 0.2 | 0.055354217756 | 338 | 21 | 239 | 78 |
| 1.2001 | 1.4001 | 0.2 | 0.015215672352 | 239 | 5 | 166 | 68 |
| 1.4001 | 1.6001 | 0.2 | 0.015510882159 | 166 | 3 | 125 | 38 |
| 1.6001 | 1.8001 | 0.2 | 0.043571462963 | 125 | 5 | 83 | 37 |
| 1.8001 | 2.0001 | 0.2 | 0.012121212184 | 83 | 1 | 61 | 21 |
| 2.0001 | 2.2001 | 0.2 | 0.015402549628 | 61 | 1 | 47 | 13 |
| 2.2001 | 2.4001 | 0.2 | 0.002976338841 | 47 | 0 | 37 | 10 |

## JSM - semi Cox


| Tstart | Thoriz | Dt | Brier | n_risk | n_event | n_alive | n_unknown |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | 
| 0.0001 | 0.2001 | 0.2 | 0.13216666217 | 4993 | 2010 | 2013 | 970 |
| 0.2001 | 0.4001 | 0.2 | 0.11480553969 | 2013 | 326 | 1162 | 525 |
| 0.4001 | 0.6001 | 0.2 | 0.08156451555 | 1162 | 118 | 744 | 300 | 
| 0.6001 | 0.8001 | 0.2 | 0.05901066681 | 744 | 56 | 488 | 200 | 
| 0.8001 | 1.0001 | 0.2 | 0.04707375326 | 488 | 23 | 338 | 127 | 
| 1.0001 | 1.2001 | 0.2 | 0.05392240669 | 338 | 21 | 239 | 78 |
| 1.2001 | 1.4001 | 0.2 | 0.01689360611 | 239 | 5 | 166 | 68 | 
| 1.4001 | 1.6001 | 0.2 | 0.01592996030 | 166 | 3 | 125 | 38 |
| 1.6001 | 1.8001 | 0.2 | 0.04602847468 | 125 | 5 | 83 | 37 |
| 1.8001 | 2.0001 | 0.2 | 0.01243010385 | 83 | 1 | 61 | 21 | 
| 2.0001 | 2.2001 | 0.2 | 0.01528940738 | 61 | 1 | 47 | 13 | 
| 2.2001 | 2.4001 | 0.2 | 0.00000000000 | 47 | 0 | 37 | 10 | 
