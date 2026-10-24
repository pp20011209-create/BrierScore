# 右設限聯合模型與動態 IPCW Brier：報告教材

本專案整理原 `cox_exp.R` 的流程，搭配逐區塊中文解說，區分「資料生成分布」與「配適模型」。僅討論右設限；不是原論文所有表格的重現。

## 從哪裡開始

- **[REPORT.md](REPORT.md)**：報告用全文；每個主要程式區塊上方都有邏輯與符號說明，並附完整程式。
- **[BASELINE_SWITCH.md](BASELINE_SWITCH.md)**：更換 exp、Weibull、lognormal、log-logistic 生成基線要改哪裡，含公式與參數表。
- **[MODEL_STATUS.md](MODEL_STATUS.md)**：原碼哪些模型實際完成，哪些名稱不符，哪些尚未完成。
- **[run_analysis.R](run_analysis.R)**：完整執行入口；會跑全部時間網格，不先做單一時間示範。
- **[R/](R/)**：依基線、生成、配適、預測與評分分成四個程式檔。

## 目前可做與不可宣稱的事

| 部分 | 狀態 |
|---|---|
| 四種 PH 基線的生成程式 | 已撰寫通用介面；保留積分＋求根；擴充分支尚未執行驗證 |
| exp 與 Weibull 配適 | 使用 `JM::jointModel()`；exp 固定 `scaleWB=1` |
| 原作者 lognormal 節點分支 | 保留為樣條 PH，不是嚴格 Lognormal-PH |
| 原作者 log-logistic 節點分支 | 本次補齊，保留為樣條 PH，不是嚴格 Loglogistic-PH |
| 半參數 Cox / JSM 分支 | 尚未納入；不能說只剩 log-logistic 一個模型 |
| 動態 Brier | 移動起點、固定跨度、條件反向 KM 的 IPCW |
| IBS、逐人權重明細、重複模擬 | 暫不執行 |
| 獨立樣本評估、數值與統計測試 | 依要求暫不加入／執行，不宣稱已驗證 |

## 執行環境

需要 R 與 `MASS`、`nlme`、`survival`、`JM`、`JMbayes`。保留原作者的 JMbayes 介面，未替換為 JMbayes2。若尚未安裝，先在 R 執行下面的安裝指令；套件版本與平台可能影響安裝方式，執行入口不會自動安裝套件。

```r
install.packages(c("MASS", "nlme", "survival", "JM", "JMbayes"))
```

將 R 工作目錄設為本資料夾，再執行下列指令。這會生成 5,000 人、配適四個模型（其中兩個為 MCMC），再跑完整時間網格，可能需要較長時間及較多記憶體。若已有對話中配適完成的物件，可閱讀 REPORT 的「接續現有工作階段」區塊，不必重新生成及配適。

```r
source("run_analysis.R", encoding = "UTF-8")
```

執行結果寫入 `results/`，預設被 `.gitignore` 排除；不要將既有論文全文、病人原始資料、套件安裝目錄或帳號憑證一併上傳。此資料夾只含教學程式與文件，沒有複製論文全文或原始臨床資料。

## GitHub 狀態

目標為使用者指定的 [BrierScore repository](https://github.com/pp20011209-create/BrierScore)，使用者指定保留私人設定。本教材不修改 repository 可見性。未自行加入 LICENSE：授權需由有權授權者決定。

## 文件依據

模型與介面依據列於 [REPORT.md](REPORT.md) 最後的套件文件區塊。原碼作為流程參考；已知不一致之處不以改名掩蓋。尚未執行本次新整理程式，不附虛構的測試通過或收斂宣告。
