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
