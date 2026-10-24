# 固定真實基線為 exp：逐段生成、配適與動態 Brier 評分

這是額外新增的獨立教學版本。原本根目錄的通式版仍保留；這裡不呼叫通式版的生成或配適函數，也不修改它們。完整可執行程式是 [analysis.R](analysis.R)，本頁則把同一份程式放在解說後面，依序讀、依序執行即可。

本版的研究問題是：**先生成一份真實基線為 exp 的資料，再把同一份資料交給不同模型，評估各自的動態 Brier score。** 不會每配適一個模型就重新抽一份資料，也不會隨機選擇生成基線。

依你的選擇，lognormal／loglogistic 保留原作者的做法：先用該分布的分位數選節點，再配適樣條 PH。因此完整比較的四個名稱是：

| 表格名稱 | 真正配適的模型 |
|---|---|
| `exp_PH` | Weibull PH 固定形狀為 1，等價於 exp 基線 |
| `weibull_PH` | Weibull PH，形狀由資料估計 |
| `spline_PH_lognormal_knots` | lognormal 分位數選節點的樣條 PH |
| `spline_PH_loglogistic_knots` | loglogistic 分位數選節點的樣條 PH |

後兩項不能寫成已完成嚴格的 Lognormal-PH／Loglogistic-PH。前兩項是 JM 的最大概似估計，後兩項是 JMbayes 的貝氏估計；比較時基線彈性和估計方法都不同，不能把差異全部歸因於分布名稱。

## 怎麼使用，才不影響原本的工作

建議開新的 R 工作階段，將工作目錄設在 BrierScore repository 根目錄。若要在目前工作階段執行，可用下列方式把本版的變數存到獨立環境 `exp_demo`，避免覆蓋你原本的 `sim`、`long`、`jointfit_exp`。`source` 會直接執行全部模型與時間點，不是先試一組。套件載入與亂數種子仍影響同一個 R 工作階段，因此新工作階段最乾淨。

```r
exp_demo <- new.env(parent = globalenv())
source("exp-only/analysis.R", local = exp_demo, encoding = "UTF-8")
exp_demo$comparison
```

若你選擇逐段複製下方程式，請在新的 R 工作階段依序貼上，不要又完整 `source` 一次造成重複配適。這次依要求不執行測試或模型運算；下列輸出名稱是程式設計，不是宣稱已取得結果。

## 01｜載入套件，定義基線風險

### 這段要做什麼

先定義「所有人共用、尚未乘上個人指標效果」的基線風險 $h_0(t)$。沿用你貼的函數名稱 `haz0_weibull`：`haz` 表示 hazard，`0` 表示 baseline，`weibull` 表示公式形式。這只是作者取的名稱，不是 R 保留字。

$$h_0(t)=\frac{k}{a}\left(\frac{t}{a}\right)^{k-1},\qquad k=\texttt{shape},\ a=\texttt{scale}.$$

本版生成資料時固定 `shape=1`，所以 $h_0(t)=1/a$；`scale=100` 對應基線風險率 0.01。這不表示每位病人的風險都不隨時間改變，因為稍後還會乘上 $\exp\{\gamma Z_i+\alpha m_i(t)\}$。

### 每個主要符號

- `haz0_weibull <- function(...)`：把一個函數存進 `haz0_weibull`。左側是函數的名字，右側是做法；此時還沒有生成病人。
- `t`：時間，可以一次傳入多個時間；`tt`：把太小的時間換成 `eps` 後的暫存結果。
- `eps=1e-12`：$10^{-12}$。`pmax(t, eps)` 逐元素取較大值，不是把整組時間只留最大值。
- `shape-1=0` 時，時間部分的零次方為 1，因此這個替換不影響本版 exp 基線。但若未來 `shape!=1`，它會改變極靠近 0 的原公式，不能當作完全精確的 Weibull 原點處理。
- `stop(...)`：發現不合法輸入就中止。這是執行時檢查，不是另外跑一份測試資料。

`MASS` 生成相關的隨機效果；`nlme` 配適指標的混合模型；`survival` 處理存活資料、Cox、參數存活模型與 KM；`JM`、`JMbayes` 配適聯合模型。本版保留 JMbayes，不改用介面不同的 JMbayes2。呼叫重要函數時加 `套件::函數`，避免同名函數被後載入的套件遮蔽。

```r
library(MASS)
library(nlme)
library(survival)
library(JM)
library(JMbayes)

haz0_weibull <- function(t, shape, scale, eps = 1e-12) {
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0 ||
      length(scale) != 1L || !is.finite(scale) || scale <= 0 ||
      length(eps) != 1L || !is.finite(eps) || eps <= 0 ||
      any(!is.finite(t)) || any(t < 0)) {
    stop("時間須非負且有限，shape、scale、eps 須為有限正數。")
  }
  tt <- pmax(t, eps)
  (shape / scale) * (tt / scale)^(shape - 1)
}
```

## 02｜cox_exp：由累積風險門檻找出事件時間

### 為什麼先積分，再找時間

個人的瞬間風險為 $h_i(u)=h_0(u)\exp\{\eta_i(u)\}$；累積風險為

$$A_i(t)=\int_0^t h_i(u)\,du,\qquad S_i(t)=\exp[-A_i(t)].$$

先抽 $U\sim U(0,1)$ 再解 $A_i(T_i)=-\log U$ 可以生成事件時間。程式等價地只抽一次 `E <- rexp(1)`，亦即 $E\sim\mathrm{Exp}(1)$；`E` 是**累積風險門檻**，不是直接抽好的事件時間。

`A` 儲存「輸入 t 就回傳累積風險」的函數，`f` 儲存 $f(t)=A(t)-E$。找 `f(t)=0` 就是在找門檻首次被累積風險碰到的時間。因為風險非負，累積風險不會隨時間下降；本版正風險下，有限解唯一。

### 先說明傳入值和內層函數

| 名稱 | 作用 |
|---|---|
| `id` | 這次在替哪位病人生成時間 |
| `eta_fun(u,id)` | 在時間 u 讀取那位病人的線性預測值 |
| `haz0_fun(u)` | 基線風險函數；未提供時使用前面的 Weibull 公式 |
| `NULL` | 尚未提供／沒有這個物件；不是數值 0，也不是無限大 |
| `A_inf` | $A_i(\infty)$ 的正確上限資訊；第三段會逐人計算並傳入 |
| `t_upper_init` | 初次搜尋右端點；函數預設 10，本版呼叫時傳入 30 |
| `t_upper_max` | 搜尋右端點最多到哪裡；不是事件時間的理論最大值 |
| `tU` | 搜尋時目前的右端點，U 代表 upper |
| `fu` | $f(tU)$，用來判斷門檻是否已被跨過 |
| `ans` | 求根得到的答案，最後回傳給外層 |

`haz0_fun <- function(u) ...` 把外層傳入的 `shape`、`scale` 記住；內層 `A` 也記住這次的 `id`、`eta_fun`。因此每次呼叫 `cox_exp`，都會建立對應那位病人的計算環境，不會將所有人混成同一條風險曲線。

### 「找不到，所以乘 2」究竟做了什麼

假設初始上限是 10。先算 $f(10)$：

- 若 $f(10)\ge0$，門檻已在 0 到 10 之間被跨過，停止擴張。
- 若 $f(10)<0$，累積風險還不夠；將 `tU` 改成 20，再算 $f(20)$。
- 若仍小於 0，繼續試 40、80……，直到跨過門檻，或碰到 `t_upper_max`。

例如假設 $E=1$、$A(10)=0.4$、$A(20)=0.8$、$A(40)=1.3$，就先找出含解的範圍，再請 `uniroot` 找到真正的時間。這些數字僅為示意，不是本程式已跑出的結果。

乘 2 是快速擴張搜尋範圍的策略，不是統計公式，也不是說病人的時間會以 2 倍成長。`min(tU*2,t_upper_max)` 防止最後一次超過指定界線。本版傳入初始上限 30，因此實際候選上限為 30、60、120……。`E` 不會在迴圈裡重抽，否則每次都換目標，便不是在求同一個人的同一個事件時間。

### 有沒有「拿上一段繼續積分」？

**你貼的版本及此版都沒有快取上一段面積。** 上限從 10 擴到 20 時，程式重新呼叫 `A(20)`，重新積分 $\int_0^{20}h_i(u)du$。沿用的是先前的搜尋上限數值 `tU`，不是先前的積分結果。最後仍把 `c(0,tU)` 交給求根器，沒有改成 `c(上一個tU,新tU)`。

數學上確實可以寫成 $A(20)=A(10)+\int_{10}^{20}h_i(u)du$，也可以利用已知符號縮小求根範圍；那是另一種加速實作。本版保留你熟悉的形式，避免一邊解釋原碼一邊換成不同演算法。不能只積分 10 到 20 的新面積，再直接拿它和整個門檻 E 比，因為那會漏掉先前的累積風險。

### integrate 的「分區間」與精度

`integrate` 會自適應細分積分區間；`subdivisions=2000L` 是允許的子區間數上限，不是預先切成 2,000 個等距時間點。`L` 表示整數。`rel.tol` 是要求的相對精度；這裡也明列 `abs.tol=rel.tol`，與原函數的預設一致。`$value` 取出積分估計值。這些是數值精度設定，不是研究的觀測時間間隔或樣本數。[R integrate 文件](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/integrate.html)

### uniroot 如何接手

抽到 $E>0$ 時，$f(0)=-E<0$；搜尋完成時 $f(tU)\ge0$。連續函數在這個範圍內有零點，`uniroot` 以內部求根演算法縮小範圍，並非固定每次只做簡單二分。`$root` 取時間答案，`tol=1e-8` 是求根的容許精度，與積分容許精度用途不同。本版補 `check.conv=TRUE`，避免求根未收斂卻只警告並繼續使用近似答案。[R uniroot 文件](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/uniroot.html)

若已知 $E\ge A_\infty$，直接回傳 `Inf` 才有模型依據。若只是在 `t_upper_max` 前沒找到，便中止，而不是把找不到解硬改成 `Inf`。`A_inf=NULL` 只表示不知道上限；`A_inf=Inf` 才表示知道累積風險沒有有限上限。積分仍可能在極端參數或很大區間出現數值問題，增加上限並不能保證解決；本次未做數值驗證。

以下保留原始主流程，只補輸入檢查、明列絕對積分容許誤差，以及求根收斂／結果檢查。

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
    haz0_fun <- function(u) haz0_weibull(u, shape, scale, eps = eps)
  }
  if (!is.function(haz0_fun) || !is.function(eta_fun)) {
    stop("haz0_fun 與 eta_fun 必須是函數。")
  }
  if (length(t_upper_init) != 1L || !is.finite(t_upper_init) || t_upper_init <= 0)
    stop("t_upper_init 必須是單一、有限且大於 0 的數值。")
  if (length(t_upper_max) != 1L || !is.finite(t_upper_max) || t_upper_max < t_upper_init)
    stop("t_upper_max 必須有限，且不得小於 t_upper_init。")
  if (length(rel.tol) != 1L || !is.finite(rel.tol) || rel.tol <= 0 ||
      length(subdivisions) != 1L || !is.finite(subdivisions) ||
      subdivisions < 1 || subdivisions != floor(subdivisions))
    stop("rel.tol 須為正數，subdivisions 須為正整數。")
  if (!is.null(A_inf)) {
    if (length(A_inf) != 1L || is.na(A_inf) || A_inf <= 0)
      stop("A_inf 必須是正數或 Inf。")
  }

  A <- function(t) {
    if (length(t) != 1L || !is.finite(t) || t < 0)
      stop("A(t) 的 t 必須是單一、有限且非負的數值。")
    if (t == 0) return(0)
    val <- integrate(
      function(u) haz0_fun(u) * exp(eta_fun(u, id)),
      lower = 0, upper = t,
      rel.tol = rel.tol, abs.tol = rel.tol,
      subdivisions = subdivisions, stop.on.error = TRUE
    )$value
    if (!is.finite(val) || val < 0)
      stop("病人 ", id, "：累積風險計算失敗。")
    val
  }

  E <- rexp(1)
  if (!is.finite(E) || E <= 0) stop("抽到無效的累積風險門檻。")
  if (!is.null(A_inf) && E >= A_inf) return(Inf)
  f <- function(t) A(t) - E
  tU <- t_upper_init
  fu <- f(tU)
  while (fu < 0) {
    if (tU >= t_upper_max) {
      stop("病人 ", id, "：搜尋上限內仍未找到事件時間；",
           "請檢查搜尋範圍、積分結果與累積風險上限。")
    }
    tU <- min(tU * 2, t_upper_max)
    fu <- f(tU)
  }
  ans <- uniroot(f, interval = c(0, tU), tol = 1e-8, check.conv = TRUE)$root
  if (!is.finite(ans) || ans <= 0)
    stop("病人 ", id, "：有限事件時間計算失敗。")
  ans
}
```

## 03｜simulate_cox_exp：每位病人的指標、事件與右設限

### 先生成個體差異，再把個體差異帶進風險

`n=5000` 是病人人數。`Z` 存每個人的不隨時間改變共變數；預設 `rbinom(n,1,0.5)` 生成 0 或 1。`b` 是 n 列、2 欄的矩陣，每列一位病人。`MASS::mvrnorm` 依 `D` 生成相關的隨機截距 `b0`、隨機斜率 `b1`。

`D` 的兩個變異數都是 0.5，標準差為 $\sqrt{0.5}\approx0.7071$；共變異數 0.01 對應相關係數 0.02。不是把 0.01 當成相關係數。`byrow=TRUE` 表示把提供的四個數字按列填進矩陣。

$$m_i(t)=(2+b_{0i})+(-0.5+b_{1i})t,\quad y_i(t)=m_i(t)+\epsilon_i(t),\quad\epsilon_i(t)\sim N(0,0.1^2).$$

`m_it` 儲存真實軌跡函數，`eta_fun` 儲存 $\eta_i(t)=\gamma Z_i+\alpha m_i(t)$。其中 `beta_fix` 是指標的固定截距／斜率，`alpha=2.5` 是指標對風險的關聯強度，`gamma=1` 是 Z 對風險的係數，三者不能混用。指標誤差不直接放進事件生成公式。

### 為什麼這裡要另外算 A_inf

本版固定 `shape=1` 且指標為直線。因此令

$$\eta_{0i}=\gamma Z_i+\alpha(\beta_0+b_{0i}),\quad s_i=\alpha(\beta_1+b_{1i}),\quad q_i=\exp(\eta_{0i})/\texttt{scale},$$

就有 $h_i(t)=q_i e^{s_i t}$。`eta0`、`slope` 分別存前兩項，且

$$A_i(\infty)=\begin{cases}q_i/(-s_i),&s_i<0,\\ \infty,&s_i\ge0.\end{cases}$$

程式用 `log_A_inf=eta0-log(scale)-log(-slope)` 再取指數，減少中間運算的數值範圍問題；若仍超出可表示範圍便中止。此公式只用來正確判斷是否有有限事件時間，有限時間仍交給前段積分＋求根。這解決先前負斜率病人可能永遠碰不到門檻的情況。

`T_true <- rep(NA_real_,n)` 先準備 n 個待填位置；`for(i in seq_len(n))` 逐人填入。`NA_real_` 是數值型的缺值，`Inf` 則是模型下無有限事件時間，兩者不一樣。有限 $A_\infty$ 時，$P(T_i=\infty\mid b_i,Z_i)=e^{-A_\infty}$，因此有些 `Inf` 可以是這個模型的合理結果。

**不要在此版只把 shape 改成 2。** 此上限公式依賴 exp 基線，所以本版會檢查 `shape==1`；這是固定 exp 情境版，換生成基線請讀原本的 [通式版切換說明](../BASELINE_SWITCH.md)。

### 右設限怎麼建立

`Cind` 是獨立抽出的隨機設限時間，`cens_rate=1.6` 是它的 rate；其平均為 $1/1.6=0.625$，不是平均 1.6，也不是最後設限比例。`admin_tau=30` 是行政截止時間。

$$C_i=\min(C_{ind,i},30),\quad X_i=\min(T_i,C_i),\quad\delta_i=I(T_i\le C_i).$$

`pmin` 是逐人取較小值，所以 `time` 是 X、`event` 是 δ。`event=0` 時只知道病人在 X 之前沒事件，並不知道真正事件何時發生，這就是右設限。它與預測起點是否等於 0 無關。

### 為什麼還要一個量測迴圈

每人的 `time` 不同，能保留的量測筆數也不同。`visit_grid` 是 0、0.2……2 共 11 個預定量測時點；`keep_rule` 決定保留 `< time[i]` 或 `<= time[i]`。最後量測時間 2 並不是行政設限時間 30。

`t_i` 是第 i 人仍可觀察的時點，`mu` 是這些時點的真實指標值，`y_i` 加上量測誤差。`long_list[[i]]` 儲存一人的資料表，`do.call(rbind,long_list)` 再按列合併。`next` 表示此人沒有可保留的量測時，跳到下一人。

回傳的是清單：`long` 多人多列、`surv` 每人一列、`censoring_rate` 為實際設限比例、`settings` 為設定。`surv$Ttrue` 和 `surv$C` 是模擬時的額外資訊，後續配適與 IPCW 不拿真實事件時間作弊。

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
    keep_rule = c("<", "<=")) {
  Z_dist <- match.arg(Z_dist)
  keep_rule <- match.arg(keep_rule)
  stopifnot(length(n) == 1L, is.finite(n), n >= 2, n == floor(n),
            length(shape) == 1L, is.finite(shape), shape == 1,
            length(scale) == 1L, is.finite(scale), scale > 0,
            length(beta_fix) == 2L, all(is.finite(beta_fix)),
            is.matrix(D), all(dim(D) == c(2, 2)), all(is.finite(D)),
            isSymmetric(D), min(eigen(D, symmetric = TRUE)$values) > 0,
            length(sigma_y) == 1L, is.finite(sigma_y), sigma_y > 0,
            length(alpha) == 1L, is.finite(alpha),
            length(gamma) == 1L, is.finite(gamma),
            length(admin_tau) == 1L, is.finite(admin_tau), admin_tau > 0,
            length(cens_rate) == 1L, is.finite(cens_rate), cens_rate > 0,
            length(visit_grid) > 0L, all(is.finite(visit_grid)),
            visit_grid[1] == 0, all(diff(visit_grid) > 0))

  Z <- if (Z_dist == "binom") rbinom(n, 1, 0.5) else rnorm(n)
  b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
  colnames(b) <- c("b0", "b1")
  m_it <- function(t, i) {
    (beta_fix[1] + b[i, "b0"]) + (beta_fix[2] + b[i, "b1"]) * t
  }
  eta_fun <- function(t, id) gamma * Z[id] + alpha * m_it(t, id)

  T_true <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    eta0 <- eta_fun(0, i)
    slope <- alpha * (beta_fix[2] + b[i, "b1"])
    if (!is.finite(eta0) || !is.finite(slope)) stop("風險參數不合法。")
    A_inf <- Inf
    if (slope < 0) {
      log_A_inf <- eta0 - log(scale) - log(-slope)
      A_inf <- exp(log_A_inf)
      if (!is.finite(A_inf) || A_inf <= 0)
        stop("病人 ", i, "：有限累積風險上限超出可表示的數值範圍。")
    }
    T_true[i] <- cox_exp(
      id = i, shape = shape, scale = scale, eta_fun = eta_fun,
      t_upper_init = admin_tau, A_inf = A_inf
    )
  }

  Cind <- rexp(n, rate = cens_rate)
  C <- pmin(Cind, admin_tau)
  time <- pmin(T_true, C)
  event <- as.integer(is.finite(T_true) & T_true <= C)
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
                       admin_tau = admin_tau, cens_rate = cens_rate, keep_rule = keep_rule))
}
```

## 04｜真正生成一次資料，固定這份資料作比較

前三段以定義函數為主，現在才呼叫 `simulate_cox_exp` 真正抽資料。`set.seed(2025)` 固定亂數起點，方便在相同程式與環境下重現；若改變抽亂數順序或演算法，不保證與先前對話的資料完全相同。

`sim` 是清單，`sim$long`、`sim$surv` 取出其中兩張表。`head` 只顯示前幾列，不代表只生成這幾人。`nrow(surv)` 應對應設定的 N=5000；`nrow(long)` 是量測總筆數，通常不同。

下面這段只跑一次。後續四個模型都使用同一個 `long` 和 `surv`，不在模型之間重新 `simulate_cox_exp`。若要比較你原本工作階段那份資料，應以保存的資料替代這個生成區塊，不能把本次新資料的分數與舊資料分數直接當成同資料比較。

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

## 05｜先建立指標模型與 Cox 初始物件

### 先對齊 id，為什麼不能省略

兩張表列數不同，但病人必須一致。`levels_id` 存同一套編號順序，`factor` 明確指定類別層級；`order` 排好病人與時間。`match(long$id,surv$id)` 回傳每筆量測對應到存活表的哪一列，存入 `ix`，用來檢查 Z 與追蹤時間一致。這不是把重複量測刪掉，而是避免模型把不同人的資料配在一起。

### lmefit 是哪一個模型

`y ~ time + Z` 表示固定效果包含截距、時間、Z；截距預設存在。`random=~time|id` 表示每人有隨機截距與隨機時間斜率。`|id` 是按病人分組，不是做「或」運算。`lmefit` 是這個初始縱向模型的估計結果。

生成端未讓 Z 直接影響指標，因此縱向 Z 係數真值為 0；仍沿用原作者模型包含 Z 的形式。`na.action=na.fail` 避免默默刪除缺值列。`opt="optim"` 指定初始混合模型的最佳化器，與事件生成時的 `uniroot` 完全不同。

### coxfit_weibull 這個名字容易誤解

`coxph(Surv(time,event)~Z)` 在這裡提供聯合模型所需的存活子模型結構與初始估計，並沒有指定 Weibull 基線。左側的名稱雖保留 `coxfit_weibull`，真正指定 Weibull 的地方在下一段 `jointModel`。

`cluster=id` 保留病人群聚設定；它影響 Cox 的穩健變異數處理，不代表後續聯合模型自動沿用同一套穩健標準誤。`x=TRUE`、`model=TRUE` 保留設計矩陣及模型資料，讓後續函數讀取。

`coxfit_exp <- coxfit_weibull` 是讓兩個聯合模型使用同一個初始 Cox 物件，不表示 exp 與 Weibull 最終配適也一樣。聯合模型會重新共同估計參數，報告要看最後的 joint 模型摘要，不只是 `lmefit` 摘要。

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

## 06｜配適 exp 和 Weibull：只差是否固定形狀

兩個呼叫都用同一份資料的 `lmefit`、Cox 初始物件，且都使用 `parameterization="value"`：風險與當下的潛在指標值 $m_i(t)$ 關聯。`timeVar="time"` 告訴套件哪一欄是縱向時間。

`method="weibull-PH-aGH"` 指定 Weibull PH 聯合模型與自適應高斯–埃米特數值積分。配適時需把每人未觀察的隨機效果積分掉；aGH 用於近似這個概似積分。它不是第二段為生成 T 而計算的時間積分，積分變數與用途都不同。

exp 那一個加上 `scaleWB=1`；在此 PH 實作，固定相關 Weibull 參數為 1 就得到常數基線。它不是生成端的 `scale=100`，也不是把 `lmefit` 的殘差標準差固定為 1。Weibull 那一個不固定，交給資料估計形狀。`jm_control` 讓兩邊使用相同控制設定。[JM jointModel 文件](https://search.r-project.org/CRAN/refmans/JM/html/jointModel.html)

`jointfit_exp`、`jointfit_weibull` 分開儲存，避免後面把同一個模型重複拿去算兩欄。檢查收斂與固定參數屬於執行時基本檢查，不等於模型已充分驗證。若形狀估計接近 1，兩個模型分數可以非常接近；單次樣本也不保證 exp 每一列都較低。

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

## 07｜原作者 lognormal 分支：選節點，再配樣條 PH

這段有兩個不同層次，不要把第一個分布名稱誤認為最後模型的基線限制。

先由 `survreg(...~1,dist="lognormal")` 對觀察存活資料配一個只有截距的參數存活模型。`ln0` 儲存結果；`mu_hat_ln` 是 log 時間的位置估計，`sd_hat_ln` 是 log 時間的尺度估計。這個初步模型沒有加入 Z 或縱向指標，其目的在於挑時間節點，不是最終的聯合模型。

`probs=seq(0.1,0.9,length.out=4)` 是四個分位機率，約為 0.1、0.3667、0.6333、0.9；不是 0.1 到 0.9 每隔 0.1 一個。`qlnorm` 將它們換成時間分位數，放進 `knots_ln`。接著只保留觀察時間範圍內的有限節點，去重排序；保留數量可能少於四個。

第二步才呼叫 `jointModelBayes`。`baseHaz="regression-splines"` 允許以樣條描述基線，`knots` 告訴它內部節點位置，`ordSpline=4` 是四階、也就是三次樣條。這種節點是樣條基底結構的接點，不是把病人的事件變成區間設限，也不是每段都固定同一風險值。[JMbayes 文件](https://search.r-project.org/CRAN/refmans/JMbayes/html/jointModelBayes.html)

`ObsTimes.knots=FALSE` 沿用原作者設定：套件自動選節點時使用事件時間而非所有觀察時間；本段已明確傳入節點，不是靠這個旗標選出 lognormal 分位數。`param="td-value"` 對應當下指標值的關聯。

`n.iter=12000` 是總迭代數，`n.burnin=2000` 捨棄前段，`n.thin=5` 對剩餘迭代抽稀。設定種子方便重現，但 MCMC 執行結束與印出摘要不等於已確認收斂。此版依要求不另外執行診斷測試。

結果明確取名 `jm_spline_lnknots`：`jm` 是 joint model，`spline` 表示樣條，`lnknots` 表示 lognormal 選節點。原作者若叫它 `jm_lognormalPH`，名稱仍不會改變實際模型。

```r
ln0 <- survival::survreg(Surv(time, event) ~ 1, data = surv, dist = "lognormal")
mu_hat_ln <- as.numeric(coef(ln0))
sd_hat_ln <- as.numeric(ln0$scale)
probs <- seq(0.1, 0.9, length.out = 4)
knots_ln <- qlnorm(probs, meanlog = mu_hat_ln, sdlog = sd_hat_ln)
knots_ln <- sort(unique(knots_ln[is.finite(knots_ln) &
                                  knots_ln > min(surv$time) & knots_ln < max(surv$time)]))
if (!length(knots_ln)) stop("沒有可用的 lognormal 內部節點。")
set.seed(2025)
jm_spline_lnknots <- JMbayes::jointModelBayes(
  lmeObject = lmefit, survObject = coxfit_weibull,
  timeVar = "time", param = "td-value", baseHaz = "regression-splines",
  control = list(knots = knots_ln, ordSpline = 4, ObsTimes.knots = FALSE,
                 n.iter = 12000, n.burnin = 2000, n.thin = 5)
)
print(summary(jm_spline_lnknots))
```

## 08｜原作者 loglogistic 分支：改節點來源，保留樣條模型

沿用同一份資料。`ll0` 是只有截距的 loglogistic 存活模型；其 log 時間可表示為 location `mu_hat_ll` 加上 scale `sd_hat_ll` 乘標準 logistic 變數，因此時間分位數是

$$Q_T(p)=\exp\{\hat\mu+\hat\sigma\,Q_{\mathrm{logistic}}(p)\}.$$

這就是 `exp(mu_hat_ll+sd_hat_ll*qlogis(probs))`。這裡的 `sd_hat_ll` 沿用容易辨認的名稱，但嚴格地說它是 logistic 的 scale，不是 log 時間的標準差；標準 logistic 的變異數為 $\pi^2/3$。不要因名稱帶 sd 就當成標準差。

得到 `knots_ll` 後，篩選節點，再配相同種類的樣條 PH，結果存 `jm_spline_llknots`。這裡沒有重新生成 loglogistic 資料，生成真相仍是 exp。它與上一段主要差在節點位置，不是兩種嚴格參數基線的比較。

```r
ll0 <- survival::survreg(Surv(time, event) ~ 1, data = surv, dist = "loglogistic")
mu_hat_ll <- as.numeric(coef(ll0))
sd_hat_ll <- as.numeric(ll0$scale)
knots_ll <- exp(mu_hat_ll + sd_hat_ll * qlogis(probs))
knots_ll <- sort(unique(knots_ll[is.finite(knots_ll) &
                                  knots_ll > min(surv$time) & knots_ll < max(surv$time)]))
if (!length(knots_ll)) stop("沒有可用的 loglogistic 內部節點。")
set.seed(2025)
jm_spline_llknots <- JMbayes::jointModelBayes(
  lmeObject = lmefit, survObject = coxfit_weibull,
  timeVar = "time", param = "td-value", baseHaz = "regression-splines",
  control = list(knots = knots_ll, ordSpline = 4, ObsTimes.knots = FALSE,
                 n.iter = 12000, n.burnin = 2000, n.thin = 5)
)
print(summary(jm_spline_llknots))
```

## 09｜預測：只讀起點以前的資料，逐人取得條件存活機率

設定起點 $s=\texttt{Tstart}$、終點 $t=\texttt{Thoriz}$，要估計的是

$$\hat p_i(s,t)=\widehat P\{T_i>t\mid T_i>s,\mathcal Y_i(s),Z_i\}.$$

意思是「已知此人在 s 時仍無事件，並知道截至 s 的指標歷史，他能存活超過 t 的機率」。例如 s=0.5、t=0.7，使用量測時點 0、0.2、0.4，不能偷讀 0.6 的量測；即使最後一次量測在 0.4，`last.time` 仍是已知存活的起點 0.5。

`ids` 是這次可評分的病人編號，下一段才選出。`di` 是起點以前的量測子表；`S_hat` 準備存每人機率，並以 id 命名。`inherits` 辨認模型來自哪個套件；JM 逐人呼叫，JMbayes 使用批次呼叫，因為兩者介面及回傳格式不同。

`survfitJM(...simulate=FALSE)` 取得代入式預測，而非完整抽樣平均。JM 在估計參數下求個人隨機效果的經驗貝氏後驗眾數，再帶入預測；不能把它描述成 aGH 直接回傳個人隨機效果。JMbayes 有不同的參數估計來源，同樣選 FALSE 不代表兩套推論方式完全相同。

`sf` 是套件回傳的預測物件，`sf$summaries` 是每人的摘要清單，`m` 是某人的時間／機率表。`extract_survprob` 按欄名 `times`、`predSurv` 與終點找結果；`row` 是符合終點的列號。`abs(... )<1e-8` 容忍浮點表示造成的極小時間差，不是容忍預測錯誤。

檢查列數唯一、機率落在 0 到 1、id 能一一對應，是避免取錯列或配錯人，不是把不理想的預測改成漂亮的結果。預測失敗時停止，不刪人後繼續計分。

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
get_S_hat <- function(jointfit, long_df, ids, Tstart, Thoriz) {
  stopifnot(length(ids) > 0L, !anyNA(ids), !anyDuplicated(ids),
            Tstart > 0, Thoriz > Tstart)
  di <- long_df[long_df$id %in% ids & long_df$time <= Tstart, , drop = FALSE]
  di <- di[order(di$id, di$time), , drop = FALSE]
  if (!all(ids %in% di$id)) stop("部分病人沒有起點以前的量測。")
  stopifnot(all(is.finite(di$time)), all(is.finite(di$y)), all(is.finite(di$Z)))
  S_hat <- setNames(rep(NA_real_, length(ids)), as.character(ids))
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

## 10｜一個起點的 IPCW Brier：三種情況與 KM 權重

### 先確定誰在這次評估中

`ids <- surv_df$id[surv_df$time>Tstart]` 只選起點還被追蹤的人，再確認此人有起點以前的量測。`R` 是這次的存活資料子表，`n_risk` 是人數 $n_s$，`n_no_history` 記錄缺少指標歷史而未納入的人。若它不為 0，評估族群限於有歷史者，不能當作所有起點存活者。

每人的理想誤差為 $[I(T_i>t)-\hat p_i]^2$，但有些人在 t 以前設限，我們不知道 $I(T_i>t)$，所以用 inverse probability of censoring weighting，簡稱 IPCW。

### 三種情況如何填進公式

令 $X_i=\min(T_i,C_i)$，$\delta_i=I(T_i\le C_i)$。在本版獨立設限下，定義條件設限存活函數 $G_s(u)=P(C_i>u\mid C_i>s)$。對起點風險集中的人：

| 程式名稱 | 情況 | 真正結局 | 權重 | 評分貢獻 |
|---|---|---:|---:|---:|
| `dead` | $s<X_i\le t,\delta_i=1$ | 0 | $1/\hat G_s(X_i-)$ | $\hat p_i^2/\hat G_s(X_i-)$ |
| `alive` | $X_i>t$ | 1 | $1/\hat G_s(t)$ | $(1-\hat p_i)^2/\hat G_s(t)$ |
| `unknown` | $s<X_i\le t,\delta_i=0$ | 未知 | 0 | 0 |

因此

$$\widehat{BS}(s,t)=\frac1{n_s}\sum_{i\in\mathcal R_s}\left[\frac{I(X_i\le t,\delta_i=1)\hat p_i(s,t)^2}{\hat G_s(X_i-)}+\frac{I(X_i>t)(1-\hat p_i(s,t))^2}{\hat G_s(t)}\right].$$

`contrib` 先全部設 0，再替前兩種人填入加權平方誤差。未知者留在分母，不能刪掉再平均；這裡分母是 `n_risk`，也不是權重總和。權重大於 1 是因為仍能觀察到的人代表了部分被設限的人。

### KM 到底估計什麼

`Surv(time,1-event)` 將設限當作 KM 事件、真正事件當作 KM 的設限，因此是反向 KM，估計設限時間的存活機率，不是疾病的存活機率。`km` 保存估計曲線；`G_step` 是依曲線查值的函數，對應數學中的 $\hat G_s$。

輸入 `R` 已限制為 $X_i>s$，所有人都活著並被追蹤到 s；在目前設限獨立的設定下，這條 KM 從條件機率 1 開始估計。因此不再額外乘一次 G(s)。如果改用全樣本 KM 估計無條件 G，才要改用 $G(u)/G(s)$ 的條件化邏輯，不可把兩種寫法混用。

`left=TRUE` 使用 `km$time<x`，取事件時間前一瞬間的 $G_s(X_i-)$；一般查終點則用 `<=`。`tail(...,1)` 取最近一個已發生的 KM 跳點，因此 KM 在兩個跳點間保持上一個值。這裡的「沿用上一段」是階梯函數查值，與前面的累積風險積分面積、搜尋上限倍增都不相同。

如果所需 G 太小，權重會很大。`min_G=0.05` 是本程式的停止門檻，不是公式中的固定常數，也不是把 G 裁切成 0.05。若評分停止，需檢查追蹤支持與時間範圍，不能默默刪除失敗者或給固定機率。

本模擬獨立抽 Cind，符合使用共同 KM 的設計。若未來設限依賴病人特徵、指標或事件過程，需要相應條件設限模型與可識別假設，不能保證同一條 KM 仍適用。暫不另輸出每人的權重明細，但計分本身已使用權重。

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

## 11｜把起點往後移：相同跨度、重新選人與讀歷史

`brier_ipcw` 做一個起點；`brier_grid` 用迴圈把所有起點做完。`Tstarts` 是起點向量，`Dt` 是固定跨度，`j` 是第幾個起點；每次令 `s=Tstarts[j]`、`t=s+Dt`。`out[[j]]` 存那次的一列表格，最後 `rbind` 合成多列表格。

例如第一次評估 0.0001 到 0.2001，下一次評估 0.2001 到 0.4001。下一次會重新選擇仍被追蹤的人，重新取得截至新起點的指標資料，重新算條件預測與 KM 權重。不是把前一次 S_hat 當成下一次答案，也不是把前一次 Brier 累加上去。

同一個 joint 模型物件一路使用，不在每個起點重新配適。因為預測的條件資訊更新，所以這仍是動態預測。`tryCatch` 在錯誤前補上起點與終點資訊，然後停止；沒有用 `NA` 跳過失敗的時間。這樣才能知道哪裡出了問題。

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

## 12｜四個模型使用相同網格，輸出可以比較的表

`fits` 是四個已配適模型的清單；`lapply` 對每個模型呼叫相同的 `brier_grid`。不是在這裡再次配適，也沒有重新生成資料。`tables` 是四張評分表；`tab_exp` 等名稱讓你能分別查看。`comparison` 再把四張表上下合併，加上 `model` 欄標明來源。

`Tstarts=seq(0.0001,2.4,by=0.2)` 產生 12 個起點，最後為 2.2001，對應最後終點 2.4001；2.4 不是一定會出現在序列中的值。0.0001 是沿用目前近零起點，並不完全等於 0。後期起點可超過最後量測時間 2，此時仍能預測，只是沒有新的指標量測可加入。

一列是「某個模型，某個起點／終點，當時合格病人的平均誤差」。`n_event+n_alive+n_unknown=n_risk`；晚期人少，分數可能起伏。Brier 越小表示在同樣評估設定下誤差越小，但不要求隨時間單調下降，真實 exp 也不保證在單次模擬的每個起點都勝過 Weibull。

輸出另存 `exp-only/results`，不覆寫原本的 `results`；`.gitignore` 排除生成的結果。`sessionInfo` 保存 R 與套件版本。存檔會覆寫此新資料夾中同名的前次結果；要保留多次模擬需另加編號或目錄。

目前為同一份資料配適與評分的樣本內表現，尚非獨立驗證或重複模擬平均。這裡也沒有計算 IBS：沿起點積分 $BS(s,s+\Delta)$ 與固定起點、沿終點積分是不同設計，需先釐清論文定義。先前縱向斜率約 -0.59 的問題仍依你的要求保留待查，不因教材整理完成就視為已解決。

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

## 報告時可以用這段串起流程

「本次固定以 exponential 基線生成一批 5,000 位病人的資料，讓個人的事件風險另外依賴 Z 與隨時間變化的潛在指標。以累積風險反解事件時間，再加入獨立右設限與量測誤差。同一份觀察資料分別配適 exp、Weibull，以及兩種選節點方式的樣條 PH 聯合模型。評估時移動預測起點、固定未來跨度 0.2，僅使用起點以前的指標歷史，最後用反向 KM 估計設限權重，計算動態 IPCW Brier score。目前是樣本內、單次模擬的比較。」

### 四種「區間／時間」不要混在一起

| 場景 | 區間／時間代表什麼 | 是否沿用前段 |
|---|---|---|
| `tU` 搜尋 | 找到含事件時間解的搜尋範圍 | 用上一個上限乘 2；未快取積分面積 |
| `integrate` | 對時間風險積分的內部自適應子區間 | 由數值積分器內部處理 |
| `knots_ln/ll` | 樣條基底的節點位置 | 不是觀察到的事件區間，也不是搜尋步長 |
| `G_step` | KM 跳點之間的階梯函數區間 | 查值時保持最近一個跳點的值 |
| `Tstart,Thoriz` | 動態預測與評分的時間窗口 | 重新選人、更新歷史；不累加前次分數 |

上述全程仍是右設限資料。預測有時間窗口，不等於資料變成區間設限。

### 資料來源與版本範圍

程式以你貼出的 `haz0_weibull`／`cox_exp` 與原 `cox_exp.R` 分支為出發點；套件介面核對來源已附在對應段落。另可參閱 [JM survfitJM](https://search.r-project.org/CRAN/refmans/JM/html/survfitJM.html)、[JMbayes 預測原始碼](https://rdrr.io/cran/JMbayes/src/R/survfitJM.JMbayes.R)、[survival coxph](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html)。本次僅撰寫與核對教材內容，不宣稱已跑過模型、數值測試或 MCMC 收斂診斷。
