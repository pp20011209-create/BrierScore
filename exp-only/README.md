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
$$h_i(t)=\text{base}_i\exp(\text{slope}_i t).$$
其中，base 是時間零的個人風險率，slope 是線性預測值隨時間變化的斜率。總累積風險為：
$$A_i(\infty)=\begin{cases}\text{base}_i/(-\text{slope}_i), & \text{slope}_i<0,\\
\infty, & \text{slope}_i\ge 0.
\end{cases}$$
將此上限傳入 cox_exp。若抽取的累積風險門檻不低於有限上限，真實事件時間記為 Inf；其餘情況使用數值積分與求根取得有限事件時間。
最後生成設限時間，建立觀察時間與事件指標：
C_i=\min(C_{ind,i},\texttt{admin_tau}),\qquad
X_i=\min(T_i,C_i),\qquad
\delta_i=I(T_i\le C_i).

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

visit_grid 為 0、0.2、…、2，共 11 個預定量測時點。keep_rule 決定是否保留結束追蹤當下的量測。t_i 是第 i 人保留下來的時間，mu 是對應的潛在指標，y_i 加入獨立常態誤差。

long_list 逐人存放量測表；每人觀察時間不同，因此列數可以不同。next 表示沒有量測可保留時跳到下一人。最後量測時間 2 與行政追蹤截止時間 30 是不同設定。

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

censoring_rate 由實際 event=0 的比例計算。settings 保存生成參數與自訂基線函數，供結果追溯。後續模型與 IPCW 使用 time、event、Z 和量測歷史，不使用 Ttrue 或 C 作為預測特徵。

```r
long_df <- if (any(!vapply(long_list, is.null, logical(1)))) {
    do.call(rbind, long_list)
  } else {
    data.frame(id = integer(0), time = numeric(0), y = numeric(0), Z = numeric(0))
  }
  surv_df <- data.frame(id = seq_len(n), time = time, event = event,
                        Z = Z, Ttrue = T_true, C = C)
  list(long = long_df, surv = surv_df, censoring_rate = mean(event == 0),
       settings = list(n = n, visit_grid = visit_grid, beta_fix = beta_fix,
                       D = D, sigma_y = sigma_y, Z_dist = Z_dist,
                       gamma = gamma, alpha = alpha, shape = shape, scale = scale,
                       admin_tau = admin_tau, cens_rate = cens_rate, keep_rule = keep_rule,
                       haz0_fun = haz0_fun))
}
```

## 04｜生成一份固定的 exp 情境資料

本段執行資料生成；前面的函數定義本身不會產生樣本。set.seed 固定亂數起點；sim 儲存回傳清單，long、surv 取出兩張表。head 只顯示前幾列；N=5000 對應 nrow(surv)，不是量測表總列數。

四個模型使用這同一份資料。更換抽亂數順序後，即使種子相同，也不應期待與舊版生成完全相同的逐人資料。本版先抽 C 再抽事件門檻，因此與舊版數值結果需分開記錄。

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

coxfit_weibull 僅是 Cox 初始物件的名稱。coxph 在這一步沒有指定 Weibull 基線；後續 jointModel 才指定基線形式。x=TRUE、model=TRUE 保留模型所需資訊；cluster=id 的 Cox 穩健變異數設定不表示最終聯合模型自動使用同樣的穩健標準誤。

coxfit_exp 使用相同初始物件，兩個最終 joint 模型仍會分別估計。

```r
lmefit <- nlme::lme(
  y ~ time + Z, random = ~ time | id, data = long,
  na.action = na.fail, control = nlme::lmeControl(opt = "optim")
)
coxfit_weibull <- survival::coxph(
  Surv(time, event) ~ Z, data = surv,
  cluster = id, x = TRUE, model = TRUE, na.action = na.fail
)
coxfit_exp <- coxfit_weibull
stopifnot(identical(unique(as.character(lmefit$groups[[1]])),
                    as.character(coxfit_weibull$model[["(cluster)"]])))
```

## 06｜exp 與 Weibull 聯合模型

parameterization="value" 表示風險連結當下潛在指標值。weibull-PH-aGH 指定 Weibull PH 及自適應高斯–埃米特積分，用於對未觀察隨機效果的概似積分；與第 02 段對時間的積分不同。

exp 配適透過 scaleWB=1 固定相關 Weibull 參數為 1，使基線為常數；這是配適模型的限制，不是資料生成演算法的解析特例，也不是生成尺度 scale=100。Weibull 配適則讓形狀由資料估計。jointfit_exp 與 jointfit_weibull 分別保存兩個結果。[JM jointModel 文件](https://search.r-project.org/CRAN/refmans/JM/html/jointModel.html)

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

ln0 是只有截距的 lognormal 存活模型，用於挑選節點。mu_hat_ln、sd_hat_ln 是 log 時間的位置與尺度。probs 是四個分位機率 0.1、約 0.3667、約 0.6333、0.9；qlnorm 換算成時間分位數。

knots_ln 僅保留觀察時間範圍內的有限節點，排序並去除重複值。此初步參數模型不含指標關聯，並非最終聯合模型。

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

baseHaz="regression-splines" 決定最終基線為樣條，knots 指定內部節點，ordSpline=4 表示三次樣條。ObsTimes.knots 是套件自動選節點時的設定；本段已明確傳入分位數節點。

n.iter 為總迭代數，n.burnin 為前段捨棄次數，n.thin 為抽稀間距。seed 在 control 明確指定套件使用的亂數種子。摘要輸出不等於完成 MCMC 收斂診斷。[JMbayes 文件](https://search.r-project.org/CRAN/refmans/JMbayes/html/jointModelBayes.html)

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

ll0 為只有截距的 loglogistic 存活模型。log 時間採 logistic 的位置–尺度表達式，時間分位數為

$$Q_T(p)=\exp\{\hat\mu+\hat\sigma\,Q_{\mathrm{logistic}}(p)\}.$$

mu_hat_ll 是位置，sd_hat_ll 雖保留 sd 名稱，實際是 logistic scale，不是標準差。qlogis 提供標準 logistic 分位數，knots_ll 儲存轉換後的時間節點。此步驟使用既有 exp 情境資料，沒有重新生成 loglogistic 資料。

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

### 09a｜依時間與欄名取預測機率

m 是某人的預測摘要表，row 找出對應終點的列號。使用 times、predSurv 欄名，避免取錯欄位；時間比較允許極小的浮點誤差。機率須有限且介於 0 與 1。

```r
# 按時間與欄名取機率，不猜「第一列第二欄」，也不刪除預測失敗者。
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

# JM 用逐人呼叫；JMbayes 保留原作者的批次預測形式。
# 本教材明確固定 simulate=FALSE：兩者都是 plug-in，但參數估計來源不同。
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

simulate=FALSE 使用代入式預測，沒有完整抽樣平均所有參數不確定性。

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

JM 分支逐人呼叫 survfitJM。此設定在估計參數下使用個人隨機效果的經驗貝氏後驗眾數進行代入式預測；配適階段的 aGH 與此處個人效果的估計步驟不同。

函數最後回傳 S_hat。預測失敗會停止，不刪除該人或填入固定機率。

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

只納入 $X_i>s$ 且起點以前有量測資料的人，形成 $\mathcal R_s$。R 是其存活資料表，n_risk 是人數 $n_s$。n_no_history 記錄因缺少歷史而未納入的人；若非零，評估族群限於具備歷史資料者。

S_hat 依 R 的 id 順序重新對齊。預測時不使用起點以後的指標，觀察結局則用於後續評分。

```r
# 一個起點、一個終點的評分：條件反向 KM + 三種 IPCW 貢獻。
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

Surv(time,1-event) 以設限作為 KM 的事件、真正事件作為 KM 的設限，估計設限存活機率。輸入已是起點後仍被追蹤者，因此在本分析的獨立設限設定下，直接估計條件 G_s，不再額外乘 G(s)。

G_step 依 KM 跳點查值：left=TRUE 使用嚴格小於，得到事件前一瞬間的左極限；一般查終點用小於等於。在相鄰跳點之間保持上一個 KM 值，這與累積風險的積分面積無關。

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

共同 KM 依賴設限獨立與權重分母為正等條件；若設限依賴病人特徵或指標歷史，需相應的條件設限模型。此分析先不另輸出逐人權重表，但計分過程仍使用上述權重。

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
# 外層迴圈沿用原作者方向：移動 Tstart，固定 Dt；不重新配適模型。
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

fits 是四個模型清單，lapply 對每個模型套用相同時間網格與評分函數。tables 保存四張評分表，comparison 上下合併並加上 model 標籤。

起點由 0.0001 開始，每次增加 0.2，共 12 個，最後起點為 2.2001，終點為 2.4001。近零起點不完全等於 0。起點超過最後量測時間 2 後仍可預測，但沒有新指標量測可加入。

每列代表一個模型在指定窗口的群體平均誤差，且 n_event+n_alive+n_unknown=n_risk。Brier 不要求隨起點單調下降；不同起點的評估族群不同，後期風險集較小。真實 exp 也不保證單次樣本的每個分數都低於 Weibull。

結果另存 exp-only/results；同名輸出在重跑時會覆寫。此處是樣本內、單次模擬評分，尚未加入獨立評估、重複模擬或 IBS。

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
rownames(comparison) <- NULL
print(comparison, digits = 10, row.names = FALSE)

output_dir <- "exp-only/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, file.path(output_dir, "brier_comparison.csv"), row.names = FALSE)
saveRDS(list(settings = sim$settings, fits = fits, tables = tables),
        file.path(output_dir, "analysis.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
```

<a id="switch-baseline"></a>
## 13｜生成基線的替換介面

本版的數值核心只需要非負、在有限追蹤區間局部可積分的基線風險，以及可在積分節點求值的 eta_fun。更換基線時保留有限追蹤、積分、求根與右設限邏輯；不需推導每個分布的無限時間累積風險上限。數值穩定性仍須依新模型檢查。傳入 haz0_fun 後，基線參數由該函數決定；simulate_cox_exp 的 shape、scale 不會自動改變自訂函數的參數。

以下是替代第 04 段的生成設定，並非在四個配適模型之間執行。生成新資料後須重新建立 long、surv 及所有配適與評分物件。函數保留 simulate_cox_exp 的歷史名稱，但可接受其他基線。

### Weibull 基線

未提供 haz0_fun 時使用 haz0_weibull。shape=1.5、scale=100 為示範設定，不代表論文設定。形狀與尺度改變後，後面的配適模型仍可保留同樣四種，比較另一個生成情境。

```r
sim <- simulate_cox_exp(n = 5000, shape = 1.5, scale = 100)
long <- sim$long
surv <- sim$surv
```

### Lognormal 基線

基線風險為密度除以存活函數。以 log 密度減 log 存活機率計算，減少尾端直接相除的數值問題；u=0 時設定其風險極限 0。meanlog=2.5、sdlog=1 是 log 時間參數。其他參數在此仍沿用本版預設；若重現不同論文情境，須一併核對 alpha、設限 rate 等設定。

```r
haz0_lognormal <- function(u, meanlog = 2.5, sdlog = 1) {
  stopifnot(all(is.finite(u)), all(u >= 0),
            length(meanlog) == 1L, is.finite(meanlog),
            length(sdlog) == 1L, is.finite(sdlog), sdlog > 0)
  out <- numeric(length(u))
  positive <- u > 0
  out[positive] <- exp(
    dlnorm(u[positive], meanlog, sdlog, log = TRUE) -
      plnorm(u[positive], meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
  )
  out
}
sim <- simulate_cox_exp(n = 5000, haz0_fun = haz0_lognormal)
long <- sim$long
surv <- sim$surv
```

### Loglogistic 基線

基線風險為 $(k/a)(u/a)^{k-1}/[1+(u/a)^k]$。正時間使用 log 尺度計算；原點依形狀指定極限。shape=2、scale=10 為示範設定。

```r
haz0_loglogistic <- function(u, shape = 2, scale = 10) {
  stopifnot(all(is.finite(u)), all(u >= 0),
            length(shape) == 1L, is.finite(shape), shape > 0,
            length(scale) == 1L, is.finite(scale), scale > 0)
  out <- numeric(length(u))
  positive <- u > 0
  z <- shape * (log(u[positive]) - log(scale))
  log_denom <- pmax(z, 0) + log1p(exp(-abs(z)))
  out[positive] <- exp(log(shape) - log(scale) +
    (shape - 1) * (log(u[positive]) - log(scale)) - log_denom)
  out[!positive] <- if (shape < 1) Inf else if (shape == 1) 1 / scale else 0
  out
}
sim <- simulate_cox_exp(n = 5000, haz0_fun = haz0_loglogistic)
long <- sim$long
surv <- sim$surv
```

上述示例只改變生成端；樣條 PH 的 lognormal／loglogistic 標籤仍是選節點方式。改變 m_it 時亦須核對縱向配適結構；有限期間積分不再依賴直線斜率公式，但配適模型不會自動跟著更改。

## 執行與範圍

從 repository 根目錄執行以下程式。獨立環境保存物件，避免覆蓋既有同名資料；套件載入與種子仍作用於同一個 R 工作階段。完整程式一次執行全部模型與時間網格。

```r
exp_demo <- new.env(parent = globalenv())
source("exp-only/analysis.R", local = exp_demo, encoding = "UTF-8")
exp_demo$comparison
```

本版未執行模型、數值測試或 MCMC 收斂診斷。有限期間方法生成的是右設限觀察資料，不提供設限者完整的潛在事件時間；需要完整 T 的後續研究須另設計生成程序。

## 方法與程式文件

- [R integrate](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/integrate.html)：自適應積分、誤差與區間設定。
- [R uniroot](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/uniroot.html)：求根區間與收斂設定。
- [JM jointModel](https://search.r-project.org/CRAN/refmans/JM/html/jointModel.html)：聯合模型方法與控制參數。
- [JM survfitJM](https://search.r-project.org/CRAN/refmans/JM/html/survfitJM.html)：條件存活預測。
- [JMbayes jointModelBayes](https://search.r-project.org/CRAN/refmans/JMbayes/html/jointModelBayes.html)：樣條基線、節點與 MCMC 設定。
