# 固定真實基線為 exp 的獨立教學版本。
# 從 repository 根目錄 source("exp-only/analysis.R", local = exp_demo)。
# 每個 BLOCK 的詳細解釋與完整程式見 exp-only/README.md。
# 此次不執行測試；結果為樣本內評分，LN/LL 分支是樣條 PH。

# ==== BLOCK 01 ====
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

# ==== BLOCK 02 ====
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

# ==== BLOCK 03 ====
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

# ==== BLOCK 04 ====
set.seed(2025)
sim <- simulate_cox_exp(n = 5000, shape = 1, scale = 100)
long <- sim$long
surv <- sim$surv
print(head(long))
print(head(surv))
print(data.frame(n_subjects = nrow(surv), n_measurements = nrow(long),
                 censoring_rate = sim$censoring_rate))

# ==== BLOCK 05 ====
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

# ==== BLOCK 06 ====
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

# ==== BLOCK 07 ====
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

# ==== BLOCK 08 ====
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

# ==== BLOCK 09 ====
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

# ==== BLOCK 10 ====
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

# ==== BLOCK 11 ====
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


# ==== BLOCK 12 ====
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
