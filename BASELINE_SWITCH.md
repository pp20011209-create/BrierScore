# 更換基線：生成資料與配適模型要分開

`baseline` 指生成端的基線分布；它不會自動改變後面採用的配適模型。研究可以設計成「同一份 exp 基線生成資料，同時配適 exp、Weibull、樣條 PH」，也可以再生成另一份 Weibull 基線資料，重做整套比較。

嚴格地說，是基線存活時間具有某種分布，並由該分布導出風險函數；不是風險函數本身是一個隨機變數且服從該分布。加入時間相依的指標後，每人的事件時間一般不再直接服從那個基線分布。

## 1. 四種基線怎麼換

令尺度為 $a>0$，形狀為 $k>0$。生成時皆使用 $h_i(t)=h_0(t)\exp\{\gamma Z_i+\alpha m_i(t)\}$。

| 生成設定 | 基線風險 $h_0(t)$，$t>0$ | 要設定的參數 |
|---|---|---|
| `exp` | $1/a$ | `scale=a`；不使用 `shape` |
| `weibull` | $(k/a)(t/a)^{k-1}$ | `shape=k`, `scale=a` |
| `lognormal` | $\phi(z)/(t\sigma[1-\Phi(z)])$，$z=(\log t-\mu)/\sigma$ | `meanlog=mu`, `sdlog=sigma` |
| `loglogistic` | $(k/a)(t/a)^{k-1}/[1+(t/a)^k]$ | `shape=k`, `scale=a` |

在 `run_analysis.R` 的 `generation` 區塊替換參數即可。以下是目前 exp 設定；`scale=100` 對應基線 rate $0.01$，不能當作 rate 100。`cens_rate` 則是設限時間的 rate，兩者用途不同。

```r
generation <- list(
  n = 5000, baseline = "exp", scale = 100,
  alpha = 2.5, gamma = 1, cens_rate = 1.6, admin_tau = 30
)
set.seed(2025)
sim <- do.call(simulate_cox_ph, generation)
long <- sim$long
surv <- sim$surv
```

若只想觀察更換成 Weibull 基線的影響，可以在生成前改以下兩個設定。`shape=1.5` 僅為示範，並非已核對的論文設定；`shape=1` 就退化成目前的 exp 基線。

```r
generation$baseline <- "weibull"
generation$shape <- 1.5
```

你貼出的 log-normal 段落使用 $\mu_0=2.5,\sigma_0=1$，指標關聯係數為 1.5，設限 rate 為 0.1。論文的指標係數叫 `beta`，目前程式叫 `alpha`；不要誤改縱向固定效果 `beta_fix`。若要依該段設定，可在生成前改如下；這會同時改變基線、關聯強度與設限強度，不能把差異全部歸因於基線。

```r
generation$baseline <- "lognormal"
generation$meanlog <- 2.5
generation$sdlog <- 1.0
generation$alpha <- 1.5
generation$cens_rate <- 0.1
```

此教材仍保留 `admin_tau=30` 的行政設限；你貼的文字只明確列出指數設限，因此不能聲稱這份設定已完整重現論文。`visit_grid` 的最大值 2 是最後量測時間，並不自動代表生存追蹤在 2 結束。

log-logistic 則改以下位置；此處形狀 2、尺度 10 只是可讀的示範值，報告若宣稱依論文模擬，必須替換成論文真正的設定。

```r
generation$baseline <- "loglogistic"
generation$shape <- 2
generation$scale <- 10
```

以上「修改設定」都必須放在 `sim <- do.call(...)` **之前**。改完要重新生成資料，再重新配適、預測、評分；已存在的 `sim` 和模型物件不會隨著清單修改而自動更新。

## 2. 配適端哪些可以沿用

| 項目 | 更換生成基線後 |
|---|---|
| 縱向資料結構、`m_it`、量測誤差 | 若研究設計沒變即可沿用 |
| 累積風險與事件時間 | 統一由 `make_baseline`、`cox_event_time` 處理 |
| `lmefit`、Cox 初始模型與 joint 模型 | 必須用新資料重配 |
| exp 配適 | `fit_jm_ph(submodels, "exp")`：Weibull PH 固定 `scaleWB=1` |
| Weibull 配適 | `fit_jm_ph(submodels, "weibull")`：不固定形狀 |
| 原作者 LN/LL 分支 | `fit_spline_ph`：只更換樣條節點的選法 |
| 嚴格 Lognormal-PH / Loglogistic-PH 配適 | 本教材尚未實作，不能用上述樣條分支冒充 |
| IPCW 公式 | 同樣的右設限假設下沿用 |
| 評估時間網格 | 檢查各時間風險集與設限存活機率是否還足夠 |

`JMbayes` 和 `JMbayes2` 是不同套件。本版保留原作者的 `JMbayes` 介面，不能只把套件名稱加上 2 就假設程式與回傳結構相同。

## 3. 為什麼不直接抽 rlnorm / rweibull

真正要反解的是 $\int_0^T h_0(u)\exp\{\gamma Z_i+\alpha m_i(u)\}\,du=E$，其中 $E\sim\mathrm{Exp}(1)$。因為 $m_i(u)$ 隨時間變化，直接抽一個基線分布的時間再乘上固定倍率，一般不符合此模型。統一積分再求根才保留同一套邏輯。

本版 `cox_event_time` 專門支援目前的直線型 $\eta_i(t)=\eta_{0i}+s_i t$。如果未來改成二次曲線或其他非線性指標，不只是替換 `m_it`：也要改事件生成函數對累積風險上限的處理。

當斜率為負，累積風險可能有有限上限。抽到的門檻超過上限，事件時間才可設為 `Inf`；搜尋到有限時間上限仍找不到根，不能直接當作永遠沒有事件。本版對 exp/Weibull 使用上限公式，對 LN/LL 使用無限區間數值積分。後者的積分誤差估計並非嚴格保證，極端參數仍可能導致數值問題；四種生成擴充本次未執行驗證。

若只需要有限追蹤期間的觀察資料，也可另設計只反解到設限時間的演算法，但那不能把未觀察到的事件時間標成真正的 `Ttrue=Inf`。本版仍保留真實事件時間欄位的語意。
