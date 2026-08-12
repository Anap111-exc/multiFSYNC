# multiFSYNC R-only 多起点验证子协议 v3 candidate2

版本日期：2026-08-12  
协议 ID：`R_ONLY_MULTISTART_V3_CANDIDATE2_20260812`  
状态：**candidate2，未冻结、未获正式执行授权**

## 1. 版本边界

本文件是 R-only 子协议的第三轮修正候选，不覆盖且不改写：

- `PROTOCOL_R_ONLY_MULTISTART_VALIDATION_V2_DRAFT_20260812.md`；
- `PROTOCOL_R_ONLY_MULTISTART_VALIDATION_V3_CANDIDATE_20260812.md`；
- v0L 已冻结结果或 R-only v1。

它必须与总协议 `PROTOCOL_V0L_V_UNSEEN_VALIDATION_V2_CANDIDATE2_20260812.md`、manifest schema `V0LV_MANIFEST_SCHEMA_V3_CANDIDATE2_20260812` 和 runner `r_route_v3_candidate2_runner.R` 同时冻结。candidate2 审阅通过只授权另建 frozen bundle，不授权生成正式数据或启动拟合。

## 2. 固定研究问题与模型

R 路径只检验真值盲优化可达性，不修改模型、先验、CAVI 或完整 ELBO。采用方案 1 的完整函数因子共享定义。正式验证使用真实 `L_f=1`、`L_s=c(1,1)`，因子内截断 `M_f=2`、`M_s=list(2,2)`；`K=5` 是 O’Sullivan 非线性惩罚块维数，不是因子数或 FPCA 截断数。

冻结路线唯一为：

`random 初始化 + 1 次公开有界预得分 + Jaoua 默认退火 c(1,1.9,100) + 普通 T=1 CAVI`。

公开预得分必须由导出的 `multiFSYNC::bayesSYNC_multi_pre_score` 执行。formal 模式禁止注入任意 `fit_function`，并验证完整源码、已安装 namespace 和公开函数签名哈希。测试/smoke 才可注入 fake。

## 3. 起点、配对与身份

每个数据固定 12 个独立 R 起点。`fit_seed`、`seed_index`、字符串 `fit_id` 和 `initialization_independence_id` 均由 candidate2 manifest 给定；同一数据内不得重复 R seed 或 independence ID。

起点 1–4 与 B0/A 使用相同核心 seed，供拟合后配对科学评价；起点 5–12 是 R 的额外独立探索起点。seed 配对从不允许复用初始化状态，各路径分别由同一整数 seed 重新初始化。

所有 terminal rows 与 `V0LV_FIT_MANIFEST_CANDIDATE2.csv` 必须按下列五个原始字段做精确字符串连接：

`fit_id + data_id + fit_seed + seed_index + initialization_independence_id`。

禁止把 `fit_id` 转成整数、模糊匹配或依赖行顺序连接。

## 4. 主端点：严格提前停止或最多 200 个普通 T=1 sweep

200 是最大普通 `T=1` 预算，不是强制所有路径跑满 200。退火 99 sweep 不计入普通 T=1 预算。

每条主路径在每个普通 sweep 后计算完整 strict practical gate。控制固定为：`min_t1=80`、短窗 20、长窗 60、连续 5 次；ELBO 绝对速率 `0.03`、每响应速率 `1e-4`、fitted NRMSE `0.001`、RSS 相对差 `0.001`、变量 PPI 最大差 `0.01`、变量 PPI 0.95 分位差 `0.01`、factor PPI 最大差 `0.01`，长窗阈值依次为 `0.003/0.006/0.01/0.01/0.01`。

- 完整 gate 首次连续通过 5 次时立即停止，状态为 `strict_practical_converged`；
- 否则运行到实际第 200 个普通 T=1 sweep，状态为 `max_budget_reached`；
- error 必须形成状态为 `error` 的完整 terminal record，不补跑、不插补。

记录实际退火 sweep、实际 T=1 sweep、首次严格收敛 sweep、终点普通 ELBO、完整短/长窗诊断和停止原因。

## 5. 严格收敛不等于 objective eligibility

两者独立记录。数值有效但到 200 尚未严格收敛的端点可以进入预注册的描述性真值盲选择；严格 gate 不是 ELBO 最大化选择的筛选条件。

R 端点 objective eligible 当且仅当：

1. 退火完整结束并进入普通 T=1；
2. 普通 ELBO 及组成项有限；
3. `objective_valid=TRUE`；
4. 没有普通 T=1 ELBO decrease/jitter 标志；
5. 输入、配置、源码和输出哈希有效；
6. terminal 状态不是 error。

若最终 winner 未严格收敛，必须写 `selected_endpoint_unfinished=TRUE`。该批可报告描述性科学结果，但不得据此称 R 路径 `production-ready`。

## 6. 真值盲 winner 与盆地支持

每个数据只在预注册 12 个 terminal rows 中选择：先保留 objective eligible 且终点 ELBO 有限的端点，再取最大普通 T=1 ELBO。完全相等时依次按较小 `fit_seed`、较小 `seed_index` 决胜。

选择不得读取或派生 true parameters、结构标签、NRMSE、ISE、载荷匹配、R/P/L、过程/得分相关、人工图形判断或辅助稳定性指标。选择表、端点表和 ELBO 盆地表先原子写入并冻结哈希，之后才可授权解封真值。

最佳 ELBO 盆地定义为距同数据最佳 eligible ELBO 不超过 `max(0.05,1e-4*abs(best_ELBO))`。报告独立 seed/independence ID 数和 endpoint-output alignment 支持；它们解释可达性，不参与 winner 选择。

## 7. warning、error、terminal 与恢复

formal 路径不得使用 `suppressWarnings()`。使用 `withCallingHandlers()` 捕获并保存每条 warning 的 condition class、message、phase、sweep、iteration、UTC 时间、elapsed time 和 call；捕获后可 muffle，但不能丢弃。

每个 terminal record 至少含：协议/runner 版本、字符串身份字段、开始/结束/运行时间、峰值内存、terminal status/reason、error class/message/call/trace 摘要、完整 warnings、实际 sweep、严格收敛、objective checks/失败原因、终点 ELBO、input/config/source/output SHA256、formal 标志和 truth-use 标志。

RDS、CSV 和 marker 均执行“同目录临时文件→反读校验→SHA256→原子 rename”；存在完整 marker 且哈希匹配时 resume 直接复用，半成品目录必须人工隔离，不能覆盖或误判成功。

## 8. 辅助科学输出稳定性

对完整 gate 以外的输出稳定性分别记录：

- `ever_reached_5_consecutive`；
- `first_reached_sweep`；
- `maximum_streak`；
- `endpoint_pass`；
- `endpoint_streak`；
- `endpoint_stable_5_consecutive`。

稳定性同时对 fitted values、RSS、变量 PPI 和 factor PPI 比较。曾经连续通过五次但终点失败必须明确显示为“曾稳定、终点不稳定”，不得压缩成一个 `stable=TRUE`。这些字段不参与 winner。

## 9. offline racing 只作真实轨迹 replay

racing 不在正式拟合中实际启用，不改变拟合、eligibility、主 winner 或成功门槛。所有主路径按第 4 节形成真实端点后再离线 replay；不要求、也不得声称全部路径都有 200 sweep 轨迹。

每条路径只使用至其实际终点的真实普通 ELBO 轨迹。screen 为 100；提前严格结束的路径作为已完成候选保留，并在后续比较中使用其真实终点 ELBO，不伪造 101–200 的诊断。其余候选按 screen ELBO 排序，至少保留 4 个，并保留距当时最佳不超过 500 的候选。

逐数据报告：主 winner 是否被保留、是否发生假淘汰、保留数、按真实终点计算的可节省 T=1 sweep。replay 的任何结果都不回写主结果。

## 10. 03/08 固定 horizon 预算审计

只对 `v0lv_03`、`v0lv_08` 的全部 12 个 R 起点审计，共 72 个阶段记录：12×2 数据×`anchor_200/to_400/to_800`。

1. 主路径若在 200 前严格停，从完整主状态复制 audit-only 分支，固定继续到累计 200；若主路径已到 200，直接登记并哈希复用为 anchor，不重复计算。
2. 从验证过的 200 状态固定继续 200 个 sweep 到累计 400，再从验证过的 400 状态固定继续 400 个 sweep 到累计 800。
3. audit 分支记录中途 gate，但不因 gate 提前停止；必须实际到达目标 horizon。
4. 主 fit error 且无可恢复完整状态时，三个审计阶段写 `audit_unavailable_due_to_main_error`，不得选择性重跑、替补或伪造 anchor。
5. 每个 horizon 先冻结真值盲 objective 排名、反事实 winner、ELBO 盆地和包含 fitted/RSS/变量 PPI/factor PPI 的 endpoint-output alignment，之后才可评价真值。
6. 200/400/800 结果都是 counterfactual audit；主 winner 和主科学结果永远来自第 4 节的原始主端点，不能被替换或重开选择。

## 11. 测试、冻结与停止门

候选必须通过解析、单元/回归测试及明确 `formal_experiment=FALSE` 的微型端到端 smoke。至少覆盖公开/历史私有预得分等价、warning/error terminal、精确字符串连接、formal fake/hash/seed 拒绝、四类输出 alignment、稳定性反例、解封前 truth lock、resume 幂等和 audit error 不重跑。

candidate2 只可在用户另行批准后复制成不可覆盖 frozen bundle；冻结与“开始正式实验”是两个独立授权。在明确正式启动前，本协议要求 `formal_generation_started=FALSE`、`formal_fit_started=FALSE`、`formal_continuation_started=FALSE`。完成 candidate2 后停止并等待审阅。

