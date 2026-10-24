# 右設限聯合模型：從生成資料到動態 Brier score
## 1. 先辨認整個研究在做什麼

流程是：設定真實模型 → 生成一批病人 → 加入右設限與指標量測 → 用觀察資料配適 → 使用起點以前的指標預測未來存活 → 用 IPCW Brier score 評估。

「生成 exp 資料」與「配適 exp 模型」是兩件事。生成端的真值由我們設定，配適端則只看觀察資料估計參數。相同的 exp 資料可以同時拿給多種模型配適；這才是同一情境下的模型比較。資料表中的 `Ttrue` 是模擬才知道的答案，本版不把它提供給模型或右設限 IPCW 評分。

目前完成的程式分支包括 exp PH、Weibull PH、兩種選節點方式的樣條 PH。因此，**原作者的 log-logistic 選節點分支可以接著做，但不能說完整研究只剩嚴格 Loglogistic-PH 配適。** 嚴格參數式 LN/LL PH 與原檔另列的半參數分支不在本版已完成範圍。

## 2. 套件、代稱與數學模型

| 套件 | 這裡的工作 |
|---|---|
| MASS | `mvrnorm` 抽每位病人的隨機截距與斜率 |
| nlme | `lme` 建立縱向混合模型，提供聯合模型初始結構 |
| survival | Cox 初始模型、`survreg` 選節點、反向 KM |
| JM | 最大概似聯合模型；本版使用 Weibull PH 與 adaptive GH |
| JMbayes | 貝氏聯合模型；本版保留原作者的樣條基線寫法 |

病人 $i$ 的潛在指標軌跡與實際量測為

$$m_i(t)=(\beta_0+b_{0i})+(\beta_1+b_{1i})t,\qquad y_i(t)=m_i(t)+\epsilon_i(t).$$

預設 `beta_fix=c(2,-0.5)`；`b` 每列是一人的兩個隨機效果，其共變異矩陣 `D` 對角線為 0.5，因此標準差是 $\sqrt{0.5}$，不是 0.5；共變異數 0.01 對應相關係數 0.02。`sigma_y=0.1` 是量測誤差標準差。

事件風險為

$$h_i(t)=h_0(t)\exp\{\eta_i(t)\},\quad\eta_i(t)=\gamma Z_i+\alpha m_i(t).$$

`gamma` 是基線共變數 `Z` 的係數，`alpha` 是潛在指標的關聯係數。風險使用沒有量測誤差的 `m_it`，不是直接使用觀察到的 `y`。

| 寫法／物件 | 讀法與用途 |
|---|---|
| `a <- ...` | 把右側結果存到左側名稱 `a`；名稱通常由作者自訂 |
| `function(...)` | 定義可重複使用的步驟，此時尚未替每人執行 |
| `shape = 1` | 呼叫函數時指定參數；與 `<-` 的用途不同 |
| `sim$surv` | 取清單 `sim` 裡名為 `surv` 的成員 |
| `x[i]`、`x[[i]]` | 前者取子集合，後者常用來取清單第 i 個內容 |
| `long[條件, ]` | 取符合條件的列，保留所有欄 |
| `&`、`&&` | 前者逐元素判斷，後者適合單一條件並短路求值 |
| `n` | 病人人數，不是縱向資料總列數 |
| `sim`、`long`、`surv` | 模擬結果清單、多人多次量測表、每人一列存活表 |
| `jointfit_*` | 已配適的模型物件，不是生成出的資料 |
| `S_hat` | 每位病人在指定起點與終點間的條件存活機率估計 |
| `tab_*` | 各起點的 Brier 評分表 |

## 3. 基線與反解事件時間

下面第一個完整區塊定義「風險如何計算」和「何時發生事件」。`make_baseline` 回傳一個清單，其中 `log_hazard` 是可在指定時間求值的函數。使用 log 風險可減少先算很小的密度再相除的數值問題；最後仍要注意指數運算的上下溢位。

`cox_event_time` 抽 $E\sim\mathrm{Exp}(1)$，等價於先抽 $U\sim U(0,1)$ 再令 $E=-\log U$。因為 $S_i(t)=\exp[-A_i(t)]$，解 $A_i(T)=E$ 就得到事件時間。`integrate` 計算累積風險，`uniroot` 找到累積風險等於門檻的時間；外層倍增搜尋只是提供求根區間。

`eta0` 是時間 0 的線性預測值，`slope` 是它隨時間的斜率；`Ainf` 是累積風險的總上限。程式區分「有理論或數值依據判定無有限事件」與「搜尋失敗」。LN/LL 的無限積分仍是數值近似，尚需後續驗證極端參數。

```r
# 四種基線：傳回 log h0(t)，以減少尾端機率相除的數值問題。
# 這裡的 scale 是時間尺度；與 JM 的 scaleWB 不是同一個參數。
make_baseline <- function(baseline = c("exp", "weibull", "lognormal", "loglogistic"),
                          shape = 1, scale = 100, meanlog = 2.5, sdlog = 1) {
  baseline <- match.arg(baseline)
  stopifnot(length(shape) == 1L, is.finite(shape), shape > 0,
            length(scale) == 1L, is.finite(scale), scale > 0,
            length(meanlog) == 1L, is.finite(meanlog),
            length(sdlog) == 1L, is.finite(sdlog), sdlog > 0)
  log_hazard <- function(t) {
    stopifnot(all(is.finite(t)), all(t >= 0))
    ans <- rep(NA_real_, length(t))
    pos <- t > 0
    u <- t[pos]
    if (baseline == "exp") {
      ans[] <- -log(scale)
    } else if (baseline == "weibull") {
      ans[pos] <- log(shape) - log(scale) + (shape - 1) * (log(u) - log(scale))
      ans[!pos] <- if (shape < 1) Inf else if (shape == 1) -log(scale) else -Inf
    } else if (baseline == "lognormal") {
      ans[pos] <- dlnorm(u, meanlog, sdlog, log = TRUE) -
        plnorm(u, meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
      ans[!pos] <- -Inf
    } else {
      z <- shape * (log(u) - log(scale))
      # log(1 + exp(z)) 的穩定寫法。
      log_denom <- pmax(z, 0) + log1p(exp(-abs(z)))
      ans[pos] <- log(shape) - log(scale) +
        (shape - 1) * (log(u) - log(scale)) - log_denom
      ans[!pos] <- if (shape < 1) Inf else if (shape == 1) -log(scale) else -Inf
    }
    ans
  }
  list(name = baseline, shape = shape, scale = scale,
       meanlog = meanlog, sdlog = sdlog, log_hazard = log_hazard)
}

# 延續原程式 integrate + uniroot 的方法，但將基線作為輸入。
# 注意：上限判斷只適用於 eta(t) = eta0 + slope*t 的直線形式。
cox_event_time <- function(eta0, slope, baseline,
                           t_upper_init = 30, t_upper_max = 1e6,
                           rel.tol = 1e-8, subdivisions = 2000L) {
  stopifnot(length(eta0) == 1L, is.finite(eta0),
            length(slope) == 1L, is.finite(slope),
            is.finite(t_upper_init), t_upper_init > 0,
            is.finite(t_upper_max), t_upper_max >= t_upper_init)
  h <- function(u) exp(baseline$log_hazard(u) + eta0 + slope * u)
  A <- function(t) {
    if (t == 0) return(0)
    val <- integrate(h, 0, t, rel.tol = rel.tol,
                     subdivisions = subdivisions)$value
    if (!is.finite(val) || val < 0) stop("有限時間的累積風險積分失敗。")
    val
  }
  E <- rexp(1)
  A_inf <- Inf
  upper_error <- 0
  if (slope < 0) {
    if (baseline$name %in% c("exp", "weibull")) {
      k <- if (baseline$name == "exp") 1 else baseline$shape
      log_A_inf <- eta0 + lgamma(k + 1) - k * log(baseline$scale) - k * log(-slope)
      A_inf <- exp(log_A_inf)
      if (!is.finite(A_inf) || A_inf <= 0)
        stop("累積風險上限超出數值範圍；請調整數值實作，不可重抽病人。")
    } else {
      # lognormal/log-logistic + 負斜率時，上限有限，改用數值積分。
      # 此誤差為積分器的估計，不是嚴格數學誤差界。
      lim <- integrate(h, 0, Inf, rel.tol = rel.tol,
                       abs.tol = 1e-10, subdivisions = subdivisions)
      A_inf <- lim$value
      upper_error <- lim$abs.error
      if (!is.finite(A_inf) || A_inf <= 0 || !is.finite(upper_error))
        stop("無法可靠計算累積風險上限。")
      if (abs(E - A_inf) <= max(10 * upper_error, 1e-10 * A_inf))
        stop("門檻太接近數值估計的累積風險上限，不能可靠判定 T 是否有限。")
    }
  }
  if (E >= A_inf) return(Inf)
  f <- function(t) A(t) - E
  tU <- t_upper_init
  fu <- f(tU)
  while (fu < 0) {
    if (tU >= t_upper_max)
      stop("搜尋上限內未找到事件時間；不能直接當作 Inf。")
    tU <- min(2 * tU, t_upper_max)
    fu <- f(tU)
  }
  ans <- uniroot(f, interval = c(0, tU), tol = 1e-8)$root
  if (!is.finite(ans) || ans <= 0) stop("事件時間計算失敗。")
  ans
}
```

## 4. 生成病人、右設限與長格式資料

下面先抽 `Z` 與隨機效果 `b`，再建立 `m_it`、`eta_fun`，逐人取得 `T_true`。原本的 `cox_exp` 角色在這版改由共用的 `cox_event_time` 處理；保留 `simulate_cox_exp` 包裝函數，方便辨認 exp 情境。

接著生成獨立設限時間 `Cind`，加入行政截止時間：$C_i=\min(C_{ind,i},\tau)$。我們觀察到 $X_i=\min(T_i,C_i)$ 與 $\delta_i=I(T_i\le C_i)$，分別存成 `surv$time` 與 `surv$event`。這就是右設限；不是區間設限，也不是將 `Tstart` 設為 0 才成為右設限。

`visit_grid` 決定原定量測時點，`keep_rule` 決定結束追蹤當下是否保留量測。`long_list` 先逐人存資料，最後合併成 `long_df`。同樣叫 `time`，`long$time` 是量測時間，`surv$time` 是事件或設限的觀察時間。`nrow(surv)` 才是你問的 N=5000；Brier 表後期的 `n_risk` 可以遠小於 5000。

```r
# 共用生成器。更換基線時改 baseline 與對應參數，不另複製四套迴圈。
# alpha 是指標對風險的係數，對應論文部分段落的 beta；不是 beta_fix。
simulate_cox_ph <- function(
    n = 5000, baseline = "exp", shape = 1, scale = 100,
    meanlog = 2.5, sdlog = 1,
    visit_grid = seq(0, by = 0.2, length.out = 11),
    beta_fix = c(2, -0.5),
    D = matrix(c(0.5, 0.01, 0.01, 0.5), 2, 2, byrow = TRUE),
    sigma_y = 0.1, Z_dist = c("binom", "norm"),
    gamma = 1, alpha = 2.5, admin_tau = 30, cens_rate = 1.6,
    keep_rule = c("<", "<="), t_upper_max = 1e6) {
  Z_dist <- match.arg(Z_dist)
  keep_rule <- match.arg(keep_rule)
  stopifnot(length(n) == 1L, is.finite(n), n >= 2, n == as.integer(n),
            length(beta_fix) == 2L, all(is.finite(beta_fix)),
            is.matrix(D), all(dim(D) == c(2, 2)), isSymmetric(D),
            all(is.finite(D)), min(eigen(D, symmetric = TRUE)$values) > 0,
            is.finite(sigma_y), sigma_y > 0,
            is.finite(gamma), is.finite(alpha),
            is.finite(admin_tau), admin_tau > 0,
            is.finite(cens_rate), cens_rate >= 0,
            length(visit_grid) >= 1L, all(is.finite(visit_grid)),
            visit_grid[1] == 0, all(diff(visit_grid) > 0))
  base_spec <- make_baseline(baseline, shape, scale, meanlog, sdlog)
  Z <- if (Z_dist == "binom") rbinom(n, 1, 0.5) else rnorm(n)
  b <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = D)
  colnames(b) <- c("b0", "b1")
  m_it <- function(t, i) {
    (beta_fix[1] + b[i, "b0"]) + (beta_fix[2] + b[i, "b1"]) * t
  }
  eta_fun <- function(t, id) gamma * Z[id] + alpha * m_it(t, id)
  T_true <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    T_true[i] <- tryCatch(
      cox_event_time(
        eta0 = eta_fun(0, i),
        slope = alpha * (beta_fix[2] + b[i, "b1"]),
        baseline = base_spec, t_upper_init = admin_tau,
        t_upper_max = t_upper_max
      ),
      error = function(e) stop("病人 ", i, "：", conditionMessage(e), call. = FALSE)
    )
  }
  Cind <- if (cens_rate == 0) rep(Inf, n) else rexp(n, rate = cens_rate)
  C <- pmin(Cind, admin_tau)
  time <- pmin(T_true, C)
  event <- as.integer(T_true <= C)
  long_list <- vector("list", n)
  for (i in seq_len(n)) {
    t_i <- if (keep_rule == "<") visit_grid[visit_grid < time[i]] else
      visit_grid[visit_grid <= time[i]]
    if (!length(t_i)) next
    mu <- m_it(t_i, i)
    y_i <- mu + rnorm(length(t_i), 0, sigma_y)
    long_list[[i]] <- data.frame(id = i, time = t_i, y = y_i, Z = Z[i])
  }
  long_df <- do.call(rbind, long_list)
  surv_df <- data.frame(id = seq_len(n), time = time, event = event,
                        Z = Z, Ttrue = T_true, C = C)
  list(long = long_df, surv = surv_df, censoring_rate = mean(event == 0),
       settings = list(n = n, baseline = baseline, shape = shape, scale = scale,
                       meanlog = meanlog, sdlog = sdlog, beta_fix = beta_fix,
                       D = D, sigma_y = sigma_y, gamma = gamma, alpha = alpha,
                       admin_tau = admin_tau, cens_rate = cens_rate))
}

# 保留熟悉的舊名稱。這個 wrapper 明確只生成 exp 基線資料。
simulate_cox_exp <- function(...) simulate_cox_ph(baseline = "exp", ...)
```

## 5. 配適：exp、Weibull 與樣條 PH

下面先對齊兩張表的病人與順序，再建立 `lmefit` 與 Cox 初始模型。`lmefit` 包含固定截距、時間、Z，以及個人隨機截距與斜率。此生成設定沒有 Z 對指標的直接效應，因此縱向 Z 係數真值為 0。最終應解讀聯合模型摘要中的縱向估計，不把初始 `lmefit` 當成最終結果。

`coxph` 本身不指定 Weibull 基線；真正開始指定 Weibull 的地方是 `jointModel(method="weibull-PH-aGH")`。此設定用自適應高斯–埃米特積分近似對隨機效果的積分。`scaleWB=1` 在這個 PH 方法中固定 Weibull 形狀為 1，形成 exp 基線；它不是生成端的尺度 `scale=100`。不固定該參數則估計 Weibull 形狀。

保留 `cluster=id` 供 Cox 物件使用，但不能把 Cox 的穩健變異數解讀成最終聯合模型也自動使用同樣的 cluster 穩健標準誤；個體內的指標相關主要由隨機效果建模。

`fit_spline_ph` 沿用原作者「先用參數分布挑節點，再配樣條」的邏輯。`survreg` 的 LN/LL 分位數只產生 `knots`；最終 `baseHaz="regression-splines"` 仍是樣條 PH。`ordSpline=4` 為四階、即三次樣條。`n.iter`、`n.burnin`、`n.thin` 控制 MCMC 迭代、暖身捨棄與抽稀；這些數值不保證收斂。此處使用 JMbayes，沒有轉成 JMbayes2。

```r
# 模型只使用觀測欄位，不把模擬真值 Ttrue、C 或 b 當成預測特徵。
prepare_submodels <- function(long, surv) {
  stopifnot(!anyNA(long$id), !anyNA(surv$id),
            !anyDuplicated(surv$id), setequal(long$id, surv$id),
            !anyDuplicated(long[c("id", "time")]),
            all(is.finite(long$time)), all(is.finite(long$y)),
            all(is.finite(long$Z)), all(is.finite(surv$Z)),
            all(is.finite(surv$time)), all(surv$time > 0),
            all(surv$event %in% c(0, 1)))
  # 同時設定相同的 factor levels，避免數字/文字排序造成配對錯誤。
  levels_id <- sort(unique(surv$id))
  long$id <- factor(long$id, levels = levels_id)
  surv$id <- factor(surv$id, levels = levels_id)
  long <- long[order(long$id, long$time), ]
  surv <- surv[order(surv$id), ]
  ix <- match(long$id, surv$id)
  stopifnot(all(long$time >= 0), all(long$time <= surv$time[ix]),
            all(long$Z == surv$Z[ix]))
  lmefit <- nlme::lme(
    y ~ time + Z, random = ~ time | id, data = long,
    na.action = na.fail, control = nlme::lmeControl(opt = "optim")
  )
  coxfit <- survival::coxph(
    survival::Surv(time, event) ~ Z, data = surv,
    cluster = id, x = TRUE, model = TRUE, na.action = na.fail
  )
  stopifnot(identical(unique(as.character(lmefit$groups[[1]])),
                      as.character(coxfit$model[["(cluster)"]])))
  list(lmefit = lmefit, coxfit = coxfit, long = long, surv = surv)
}

# exp/Weibull 都留在 JM 最大概似框架，只有 shape 是否固定的差別。
fit_jm_ph <- function(submodels, baseline = c("exp", "weibull"), control = list()) {
  baseline <- match.arg(baseline)
  fit <- JM::jointModel(
    lmeObject = submodels$lmefit, survObject = submodels$coxfit,
    timeVar = "time", parameterization = "value",
    method = "weibull-PH-aGH",
    scaleWB = if (baseline == "exp") 1 else NULL,
    control = control
  )
  if (is.null(fit$convergence) || fit$convergence != 0)
    stop("聯合模型未回報成功收斂。")
  if (baseline == "exp" && !isTRUE(all.equal(as.numeric(fit$coefficients$sigma.t), 1)))
    stop("exp 的形狀沒有正確固定為 1。")
  fit
}

# 保留原作者用分位數選節點的做法；名稱明確標成 spline PH。
# 此函數不會配適嚴格的 Lognormal-PH 或 Loglogistic-PH。
fit_spline_ph <- function(submodels, knot_dist = c("lognormal", "loglogistic"),
                          seed = 2025, n.iter = 12000, n.burnin = 2000,
                          n.thin = 5) {
  knot_dist <- match.arg(knot_dist)
  surv <- submodels$surv
  ln0 <- survival::survreg(survival::Surv(time, event) ~ 1,
                           data = surv, dist = knot_dist)
  mu_hat <- as.numeric(coef(ln0))
  sd_hat <- as.numeric(ln0$scale)
  stopifnot(is.finite(mu_hat), is.finite(sd_hat), sd_hat > 0)
  probs <- seq(0.1, 0.9, length.out = 4)
  knots <- if (knot_dist == "lognormal") {
    qlnorm(probs, meanlog = mu_hat, sdlog = sd_hat)
  } else {
    exp(mu_hat + sd_hat * qlogis(probs))
  }
  knots <- sort(unique(knots[is.finite(knots) &
                              knots > min(surv$time) & knots < max(surv$time)]))
  if (!length(knots)) stop("沒有可用的內部節點，請檢查時間範圍。")
  fit <- JMbayes::jointModelBayes(
    lmeObject = submodels$lmefit, survObject = submodels$coxfit,
    timeVar = "time", param = "td-value", baseHaz = "regression-splines",
    control = list(knots = knots, ordSpline = 4, n.iter = n.iter,
                   n.burnin = n.burnin, n.thin = n.thin, seed = seed)
  )
  fit
}
```

## 6. 動態預測與 IPCW 三種情況

固定起點 $s=\texttt{Tstart}$，終點 $t=s+\texttt{Dt}$。只納入 $X_i>s$ 且起點以前有量測的人，形成 $\mathcal R_s$，人數為 $n_s$。對每人預測

$$\hat p_i(s,t)=\widehat P(T_i>t\mid T_i>s,\mathcal Y_i(s),Z_i).$$

`get_S_hat` 只提供截至起點的歷史資料，不能偷用未來量測。模型在不同起點不重配；個人歷史與條件存活時間改變，因此預測會更新。`simulate=FALSE` 是代入式預測：JM 使用估計參數與個人隨機效果的經驗貝氏後驗眾數；JMbayes 使用其參數後驗摘要與個人隨機效果估計。不能把這一步說成再次以 aGH 直接算出每人的隨機效果，也不能說已完整平均所有參數與隨機效果的不確定性。

理想平方誤差是 $[I(T_i>t)-\hat p_i(s,t)]^2$。但有人在終點前設限，結局未知，所以需要 IPCW。令

$$G_s(u)=P(C_i>u\mid C_i>s),\qquad u\ge s.$$

在本版獨立設限設定下，以當次風險集做反向 KM：`Surv(time, 1-event)` 把設限當作 KM 的事件，把真正事件當作 KM 的設限。由於輸入已經是起點後仍被追蹤者，這裡直接估計條件 $G_s$，不再額外乘 $G(s)$。若用全體樣本估計無條件 $G$，則理論上需用 $G(u)/G(s)$ 轉換；兩種有限樣本估計不必完全相同。

| 情況 | 可知的結局 $I(T_i>t)$ | 權重 | 對評分的貢獻 |
|---|---:|---:|---:|
| $s<X_i\le t,\delta_i=1$：期間發生事件 | 0 | $1/\hat G_s(X_i-)$ | $\hat p_i^2/\hat G_s(X_i-)$ |
| $X_i>t$：追蹤超過終點 | 1 | $1/\hat G_s(t)$ | $(1-\hat p_i)^2/\hat G_s(t)$ |
| $s<X_i\le t,\delta_i=0$：期間設限 | 未知 | 0 | 0 |

因此

$$\widehat{BS}(s,t)=\frac1{n_s}\sum_{i\in\mathcal R_s}\left\{\frac{I(X_i\le t,\delta_i=1)\hat p_i(s,t)^2}{\hat G_s(X_i-)}+\frac{I(X_i>t)[1-\hat p_i(s,t)]^2}{\hat G_s(t)}\right\}.$$

`X_i-` 指事件發生前的左極限。未知者雖貢獻 0，仍留在分母 `n_risk`；不能只平均已知者，也不是除以權重總和。假設是設限足夠獨立且所需的 G 為正；若設限依賴病人特徵或歷史，通常需要相應的條件設限模型，不能直接沿用一條共同 KM。

下面 `dead`、`alive`、`unknown` 就是這三種情況。`min_G=0.05` 是本實作的停止門檻，不是 IPCW 定義，不會把 G 強制改成 0.05。`extract_survprob` 的檢查是確認抓到正確終點與合法機率，並沒有重新配適或「修正」估計。暫不輸出逐人權重表，但計算 Brier 時仍必須使用這些權重。

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

## 7. 一次跑完整情境與所有時間點

下面是分析入口的完整內容。先載入四個函數檔，再生成一次資料，配適四個明確標示的模型，最後在相同時間網格評分。若你已經有正在使用的 `sim`，應跳過入口中的生成區塊，從 `long <- sim$long`、`surv <- sim$surv` 接續，避免不小心換一份資料而拿新分數與舊分數比較。

`fits` 是模型清單，`lapply` 對清單中每個模型跑相同步驟。`Tstarts` 移動、`Dt=0.2` 固定；這裡有 12 個起點，從 0.0001 到 2.2001，最後終點為 2.4001。使用 0.0001 沿用目前流程，並不與起點 0 完全相同。超過最後量測時間 2 仍可預測，但沒有新的指標資訊可加入。

程式直接執行完整流程，不安排單一時間點試跑。N=5000、兩個 MCMC 模型與多起點預測可能需較長時間。本次只準備程式，沒有執行；模型結果與 MCMC 收斂尚待實際分析確認。`results` 保存資料設定、模型、分數與套件版本，預設由 Git 忽略。

```r
# 在 github-right-censoring 資料夾作為 R 工作目錄時執行。
# 這是完整分析入口，不是測試程式；會配適四個模型並跑完整時間網格。
# 若目前只要閱讀，請開啟 REPORT.md，不必 source 本檔。
needed <- c("MASS", "nlme", "survival", "JM", "JMbayes")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("請先安裝套件：", paste(missing, collapse = ", "))

source("R/01_baseline.R", encoding = "UTF-8")
source("R/02_simulation.R", encoding = "UTF-8")
source("R/03_fitting.R", encoding = "UTF-8")
source("R/04_prediction_brier.R", encoding = "UTF-8")

# 區塊一：生成一次資料。換生成分布的位置就在這裡。
generation <- list(
  n = 5000, baseline = "exp", shape = 1, scale = 100,
  meanlog = 2.5, sdlog = 1, alpha = 2.5, gamma = 1,
  cens_rate = 1.6, admin_tau = 30
)
set.seed(2025)
sim <- do.call(simulate_cox_ph, generation)
long <- sim$long
surv <- sim$surv
print(data.frame(n_subjects = nrow(surv), n_measurements = nrow(long),
                 censoring_rate = sim$censoring_rate))

# 區塊二：同一份資料配適四個模型，保留不同物件，不覆蓋前一模型。
submodels <- prepare_submodels(long, surv)
long <- submodels$long
surv <- submodels$surv
lmefit <- submodels$lmefit
coxfit_weibull <- submodels$coxfit
coxfit_exp <- submodels$coxfit

# 兩者採相同積分設定；若做數值精度敏感度分析，要一起明確指定。
jm_control <- list()
jointfit_exp <- fit_jm_ph(submodels, "exp", control = jm_control)
jointfit_weibull <- fit_jm_ph(submodels, "weibull", control = jm_control)
jm_spline_lnknots <- fit_spline_ph(submodels, "lognormal", seed = 2025)
jm_spline_llknots <- fit_spline_ph(submodels, "loglogistic", seed = 2025)

print(summary(jointfit_exp))
print(summary(jointfit_weibull))
print(summary(jm_spline_lnknots))
print(summary(jm_spline_llknots))

# 區塊三：每個模型都用相同的起點、跨度、病人選取與 IPCW 規則。
# 一次執行全部時間，不先安排單一時間點測試。
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
print(comparison, digits = 8, row.names = FALSE)

# 區塊四：保存可追溯的設定、模型、分數與套件版本；產物預設不放 Git。
dir.create("results", showWarnings = FALSE)
saveRDS(list(settings = sim$settings, fits = fits, tables = tables),
        "results/analysis.rds")
write.csv(comparison, "results/brier_comparison.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "results/sessionInfo.txt")

# 目前是樣本內評分；尚未完成獨立測試、重複模擬或 MCMC 收斂確認。
# 本檔不會把樣條 PH 標成參數式 Lognormal-PH/Loglogistic-PH，也不計算 IBS。
```

### 接續現有工作階段

如果目前 R 裡已有 `long`、`surv`、`jointfit_exp`、`jointfit_weibull`，不必再生成或重配這兩個模型。下面只載入新的評分函數，直接評估全部時間。樣條模型也可同樣傳入 `brier_grid`；請使用你工作階段中真正存在的物件名稱，不要把模型換個名稱就視為另一種分布。

```r
source("R/04_prediction_brier.R", encoding = "UTF-8")
Tstarts <- seq(0.0001, 2.4, by = 0.2)
Dt <- 0.2
tab_exp <- brier_grid(jointfit_exp, long, surv, Tstarts, Dt)
tab_weibull <- brier_grid(jointfit_weibull, long, surv, Tstarts, Dt)
```

## 8. 結果要如何報告

一列表示一個模型在一組起點／終點上的群體平均誤差，不是單一病人的分數。`n_event+n_alive+n_unknown=n_risk`；`n_no_history` 是因缺少歷史而未納入的起點存活者。若該數字不為 0，評估族群已限於有歷史資料的人，不能忽略這個選取條件。

Brier 越低表示該族群與時間窗口內平方預測誤差越小，但不要求隨起點單調下降。後期事件少、風險集小，估計較不穩定；不同起點比較的也不是完全相同的人。資料由 exp 基線生成，也不能保證單次模擬中 exp 的每個分數必定小於 Weibull。

本版使用同一份資料配適與評分，是樣本內表現，不能直接稱為驗證後的泛化能力。測試／獨立評估依本次要求先不做。單次 N=5000 也不是論文多次模擬後的平均與標準差。

目前還沒有算 IBS。若要對移動起點曲線 $BS(s,s+\Delta)$ 積分，必須另說明積分變數是起點 s、固定預測跨度以及積分範圍；它與固定 s、改變未來終點 t 的整合評分不是相同設計，不能直接把這 12 列平均就宣稱重現論文 IBS。

你先前提到縱向斜率約 -0.5912，依你要求暫不在這版追查；保留為尚未解決項目，不能因流程整理完成就視為已證實無誤。

## 9. 可核對的套件文件

- [JM jointModel](https://search.r-project.org/CRAN/refmans/JM/html/jointModel.html)：方法與 scaleWB 設定。
- [JM survfitJM](https://search.r-project.org/CRAN/refmans/JM/html/survfitJM.html)：條件存活預測與 simulate 選項。
- [JMbayes jointModelBayes](https://search.r-project.org/CRAN/refmans/JMbayes/html/jointModelBayes.html)：樣條基線與控制參數。
- [JMbayes survfitJM 原始碼](https://rdrr.io/cran/JMbayes/src/R/survfitJM.JMbayes.R)：預測摘要、病人命名與代入式預測。
- [survival Surv](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/Surv.html)：右設限資料格式。
- [survival coxph](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html)：Cox 模型與 cluster。

上述為套件來源；IPCW 區塊明確列出本版採用的估計式與假設，並非聲稱已逐式核對尚未完整提供的論文。
