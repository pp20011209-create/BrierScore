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
