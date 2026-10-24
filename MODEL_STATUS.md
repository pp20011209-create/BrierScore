# 原碼與本教材的模型對照

原始流程參考使用者提供的 `cox_exp.R`。本文件不包含原論文全文，以下是對程式呼叫的核對結論。

| 原碼片段／名稱 | 實際作用 | 本教材處理 |
|---|---|---|
| `haz0_weibull` + `shape=1` | exp 基線生成 | 以通用 `make_baseline("exp")` 表示 |
| `cox_exp` | 生成單人事件時間，不是配適 | 通用核心改名 `cox_event_time`，避免與配適混淆 |
| `simulate_cox_exp` | 生成一批 exp 基線資料 | 保留相容入口；底層呼叫 `simulate_cox_ph` |
| `coxfit_weibull <- coxph(...)` | Cox 初步子模型，未指定 Weibull | 說明名稱不是模型設定；保留 cluster |
| `method="weibull-PH-aGH"` | Weibull-PH 聯合配適 | 保留 JM，形狀自由估計 |
| JM `scaleWB=1` | 在此 Weibull-PH 實作中把 shape 固定為 1 | 用於 exp-PH 配適；基線大小仍估計 |
| `jointModelBayes(baseHaz="regression-splines")` + lognormal 分位數節點 | 樣條 PH，節點由初步 lognormal 配適協助決定 | `jm_spline_lnknots` |
| 同上但使用 `exp(mu + sigma*qlogis(p))` 節點 | 樣條 PH，節點由初步 log-logistic 配適協助決定 | `jm_spline_llknots` |
| `ObsTimes.knots=FALSE` | 未手動提供節點時，用事件時間決定節點 | 移除「防止覆蓋節點」的錯誤註解 |
| `data_long` / `Surv(start,stop,event)` / `JSM::jmodelTM` | 另一個多列追蹤與聯合模型分支 | 尚未移植；切區間不等於區間設限 |
| 直接取預測表第一列第二欄 | 可能拿錯時間或欄位 | 按欄名與時間擷取；按 ID 配對 |
| 預測失敗後刪人評分 | 改變評估族群，可能偏低 | 停止並報出起點，不靜默刪除 |
| 搜尋超時直接回傳 Inf | 把數值困難誤當無限事件時間 | 用模型上限判斷；無法確認就停止 |
| `compute_IBS(tab)` | 對移動起點、固定跨度的曲線做梯形平均 | 暫不執行；不可與固定起點的 IBS 混稱 |

## 不該在報告中說的話

- 「已完成四種嚴格參數式基線的配適比較」：目前只有 exp、Weibull，加上兩種選節點方式的樣條 PH。
- 「換到 JMbayes2 就是 lognormal 基線」：套件名稱不能代替基線設定。
- 「Brier 應該隨時間單調下降」：評估族群、事件與設限權重都會變。
- 「exp 真實生成，所以 exp 在每一次、每個起點一定最低」：有限樣本、估計誤差與評分目標不同，不保證該排序。
- 「5,000 是模型名稱」：它是生成時的病人人數；`n_risk` 是當次評估人數。

## 尚未解決的限制

使用者先前得到的聯合縱向時間係數約 -0.5912，與當時設定 -0.5 有明顯差距；使用者選擇暫緩處理。本教材沒有宣稱已修正此差異。MCMC 收斂、積分精度、參數回復、獨立測試資料評估與重複模擬亦未完成。

原論文敘述的重複次數、設限率與這份原碼設定不完全相同。不能拿本教材執行一次的輸出直接當作論文表格重現。
