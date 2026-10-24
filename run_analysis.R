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
