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
