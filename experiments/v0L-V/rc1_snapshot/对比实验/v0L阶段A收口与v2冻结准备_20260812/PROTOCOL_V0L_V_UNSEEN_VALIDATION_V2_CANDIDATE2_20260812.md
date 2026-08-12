# multiFSYNC v0L-V 未见数据验证总协议 v2 candidate2

版本日期：2026-08-12  
协议 ID：`V0LV_UNSEEN_VALIDATION_V2_CANDIDATE2_20260812`  
对应 R-only ID：`R_ONLY_MULTISTART_V3_CANDIDATE2_20260812`  
对应 manifest schema：`V0LV_MANIFEST_SCHEMA_V3_CANDIDATE2_20260812`  
状态：**candidate2，未冻结、未生成正式数据、未启动正式拟合或 continuation**

## 1. 独立版本线与权限

本文件是 v0L-V 总协议 v2 的第三轮候选，不把 manifest v2 当成总协议 v2，也不覆盖：

- `PROTOCOL_V0L_V_UNSEEN_VALIDATION_V1_DRAFT_20260812.md`；
- `PROTOCOL_V0L_V_UNSEEN_VALIDATION_V2_CANDIDATE_20260812.md`；
- R-only v2/v3 candidate；
- `protocol_manifests_20260812_v2/` 或 `protocol_manifests_20260812_v3_candidate/`。

本轮只交付候选协议、runner/controller/evaluator、manifest、测试和非正式 smoke。建立 frozen bundle 与开始正式实验是两个独立的后续用户授权；任何一个都不在本候选范围内。

## 2. 目标与不变量

研究目标是以新数据检验 R 路径的真值盲优化可达性，并与 B0、A 和 pooled bayesSYNC 对照。未经新证据，不改变 multiFSYNC 模型、先验、CAVI 或完整 ELBO。采用方案 1 完整函数共享；真实/拟合 `L_f=1,L_s=c(1,1),M_f=2,M_s=list(2,2)`，`K=5`。

所有拟合、停止和 winner 选择只读取 observation-only bundle；真值先独立密封。任何 true parameters、结构、NRMSE、ISE、载荷方向匹配或人工观察都不得进入拟合、停止或选择。

## 3. 十个新数据与生成契约

正式数据 ID 和 seed 固定为：

| data_id | seed | 预算审计 |
|---|---:|---|
| v0lv_01 | 82612001 | 否 |
| v0lv_02 | 82612002 | 否 |
| v0lv_03 | 82612003 | 是 |
| v0lv_04 | 82612004 | 否 |
| v0lv_05 | 82612005 | 否 |
| v0lv_06 | 82612006 | 否 |
| v0lv_07 | 82612007 | 否 |
| v0lv_08 | 82612008 | 是 |
| v0lv_09 | 82612009 | 否 |
| v0lv_10 | 82612010 | 否 |

生成器为 `simulate_multi_study_structured`：`S=2,n_s=c(30,30),p=500,d=0,K=5,n_obs=8,common_grid=FALSE,n_obs_range=c(6,9),sigma_eps=0.3,bool_sparse_loadings=TRUE,prop_sparse=0.9,score_var_decay=TRUE,n_dense=201,mean_amp=0.6,beta_amp=0,bs_degree=c(2,3),identified_loadings=TRUE,sparsity_mode="fixed"`，并使用冻结的 v0F constant-per-active-loading 校准。生成后不得更改数据分布。

密封前逐数据强制断言：

- `S,n_s,p,d,L_f,L_s,M_f,M_s,K` 和每人 6–9 个观察点精确一致；
- A、B1、B2 各 50 个活跃变量；`supp(A)` 分别与 `supp(B1)`、`supp(B2)` 不交；
- B1/B2 在同一特异变量组内各自抽样，允许支持重叠；只要求其方向绝对余弦不超过 0.8，不宣称三方向两两支持不交；
- 载荷范数、Gram、特征函数离散 Gram 正交归一、FPCA 特征值/分量方差、误差标准差和 `mean_amp=0.6` 均通过；
- dense 与 observation reconstruction 恒等式通过；
- observation-only bundle 不含 `true_params`、结构标签、生成过程/得分、评价指标或 truth hashes。

`mean_amp` 等 provenance 只进入密封 truth/provenance bundle。控制器先写 observation bundle 和 sealed truth，再写带两者 SHA256 的 data-seal marker。

## 4. 210 条主拟合、seed 配对与路径

每个数据 21 条：B0 4、A 4、R 12、pooled bayesSYNC 1；共 B0=40、A=40、R=120、pooled=10，总计 210。

- **B0**：random，无预得分、无退火、普通 T=1 最多 200；
- **A**：random，无预得分、`c(1,1.9,100)` 退火、普通 T=1 最多 200；
- **R**：random、公开接口预得分 1 次、同一默认退火、普通 T=1 最多 200；详细契约以 R-only candidate2 为准；
- **pooled bayesSYNC**：两个研究的全部受试者按 manifest 顺序合并后只拟合一次，不传 study label，`Q=3,M=2,K=5,n_g=51,anneal=c(1,1.9,100),maxit=299`。禁止分研究分别拟合。

B0/A/R 的 seed index 1–4 在同一数据内成对，整数为 `82620000+100*data_index+seed_index`；R 的 index 5–12 为额外起点。pooled seed 为 `82630000+data_index`。相同整数只表示 RNG 配对，不共享状态。所有 fit identity 必须与 manifest 精确连接。

## 5. 停止、terminal 与 objective

B0/A/R 的普通 T=1 预算均是**最多 200**。完整 strict practical gate 一旦满足立即停止；未满足才到 200 并记 `max_budget_reached`，不得强制所有路径跑满 200。退火 sweep 单独计数。strict convergence 和 objective eligibility 是独立字段；数值有效的未严格端点可按预注册规则参与描述性选择，选中时记 `selected_endpoint_unfinished=TRUE`。

pooled 只使用其参考实现的预注册停止/最大迭代规则；其 ELBO 只判断该 pooled endpoint 数值有效，不与 multiFSYNC ELBO 比较。

全部 210 个 fit 均须有原子化 terminal record；warning 使用 condition handler 完整捕获，不使用 `suppressWarnings()`。error、警告、实际 sweep、时间、峰值内存、eligibility、所有输入/配置/完整源码/输出哈希均保留。成功 marker 和 RDS 哈希同时存在才算可恢复完成；resume 不重复计算或覆盖。

## 6. 真值盲终点选择与解封顺序

每个数据分别对 B0 的 4 个端点、A 的 4 个端点和 R 的 12 个端点按“eligible 最大同模型普通 T=1 ELBO”选择；完全相等依次按较小 fit seed、较小 seed index。pooled 是单一预注册 endpoint。所有四方法选择合并表原子写入并冻结 SHA256 后，才允许用户授权解封该数据真值。

选择前严禁读取 sealed truth。正式解封同时要求：完整 all-method selection marker、selection hash、sealed truth hash和显式 operator authorization。Racing、辅助稳定、结构或科学指标不参与选择。

## 7. R 的 offline racing 与 03/08 审计

offline racing 只使用每条主路径至真实实际终点的轨迹，screen=100、至少保留 4、ELBO margin=500。提前严格结束者作为完成候选保留；不伪造其后 sweep。报告主 winner 保留、假淘汰和估算节省，不改变主拟合、winner 或成功门。

只对数据 03/08 的全部 12 个 R 起点建立 200→400→800 audit-only 固定 horizon，共 72 条阶段记录。早停主状态复制并补到累计 200；已达 200 则哈希复用；再顺序固定到 400、800。中途 strict gate 只记录、不提前停。主 error 无状态时写 `audit_unavailable_due_to_main_error`，禁止重跑或伪造。

每个 horizon 先冻结真值盲排名、反事实 winner、ELBO 盆地和 fitted/RSS/变量 PPI/factor PPI alignment，之后才可解封评价。主 winner 和主科学结论始终来自原始严格停止或最多 200 的主端点。

## 8. multiFSYNC 身份匹配与评价顺序

冻结 v0H/v0J 评价器先按载荷方向和参数块一对一匹配共享/特异身份，再固定符号并评价特征函数、因子过程、FPCA 得分和完整贡献；不得因某一后续指标重新匹配。

报告顺序固定为：

1. 观测点和 dense-grid 重构；
2. 特征函数 total ISE、O’Sullivan 投影下界、excess ISE、函数子空间及协方差核；
3. 载荷恢复；
4. FPCA 得分相关；
5. 因子过程误差和完整贡献；
6. 完整 `1/1/1`、missing/misplaced/duplicate/extra/split 与共享/特异结构；
7. 运行时间、内存、warnings、objective、严格收敛和盆地支持。

factor PPI 的唯一规则为 `>=0.5`。真实 L/M 下的保留数只作辅助结果，不宣称已完成维数选择性能比较。

## 9. pooled bayesSYNC 评价公式

令三个真载荷方向为 `u_1=A0,u_2=B10,u_3=B20`，pooled 估计方向为 `b_q,q=1,2,3`。使用全部三个候选，不按 factor PPI 筛除。匹配排列为

`pi* = argmax_pi sum_r |u_r' b_pi(r)|/(||u_r|| ||b_pi(r)||)`，

符号 `s_r=sign(u_r'b_pi*(r))`。factor PPI `>=0.5` 仅形成 `auxiliary_factor_retained`，不改变匹配或主评价。

共享方向的 active studies 是 `{1,2}`；B10 仅 `{1}`；B20 仅 `{2}`。过程相关/NRMSE、FPCA 得分相关及完整贡献 NRMSE 只在真实 active studies 连接计算。对方向 r 的完整贡献，真值为 `u_r h_r(t)`，估计为 `b_pi*(r) hhat_pi*(r)(t)`；载荷/函数符号联动对齐。

对 B10/B20，inactive-study leakage 单列为：

- `RMS_inactive = sqrt(mean(Chat_inactive^2))`；
- `energy_fraction = sum(Chat_inactive^2)/(sum(Chat_active^2)+sum(Chat_inactive^2))`。

共享方向没有 inactive study，上述字段为真正 `NA`，并写 `shared_direction_is_active_in_both_studies`。pooled 模型没有共享/特异参数块，所以 misplaced 和 block purity 为真正 `NA` 并写原因，绝不编码为 0、成功或失败。

每个匹配方向的两组 FPCA 函数在冻结 grid 上一对一匹配。报告 total sign-aligned ISE、对精确 pooled O’Sullivan 设计空间的 projection-floor ISE、`max(0,total-floor)` excess ISE，以及只在 active subjects 上的 FPCA score correlation。pooled observed/dense reconstruction 对全部受试者一次汇总。pooled 与 multiFSYNC 的 ELBO 绝对值严禁跨模型比较。

## 10. 统计单位和成功门

保留全部逐起点原始结果。核心四 seed 的 A−B0、R−A、R−B0 在同数据同 seed 配对；误差类定义为 candidate−comparator，相关/recall/purity 类反向定义，使差值 `<0` 始终表示 candidate 更好。每数据先取四对的中位差，再把 10 个 data_id 作为独立 cluster。

报告原始分布、中位数、IQR、去除零差后的双侧精确符号检验和 data-cluster bootstrap percentile 95% 区间；`B=9999`、seed `82640001`，同一 bootstrap index 保持预注册对比配对。

跨数据中位差 `<=0` 只称**中位无恶化/优势门槛**，不是带 margin 的非劣效检验。恶化数据数保留为描述性计数，不是独立成功门，不发明非劣界值。

预注册门：

1. 每数据 12 个 R terminal、至少一个 eligible，得到 10 个真值盲 R winner；
2. 解封后至少 8/10 R winner 完整正确 `1/1/1`；
3. R−A、R−B0 的 `missing+misplaced` 跨数据中位差 `<=0`；
4. 观测 NRMSE、dense NRMSE、feature excess ISE 的相同对比逐项中位差 `<=0`；
5. 03/08 共 72 个审计 terminal 到齐或明确记录主 error 不可用，且不替换主端点；
6. 盆地与独立 alignment 支持完整披露；不足时 R 保持 diagnostic，不追加阈值；
7. 任一所选 R endpoint unfinished 时，本批不得称 production-ready；
8. 10 个 pooled terminal 全部登记，error/warning 不插补。

## 11. 完整源码绑定与正式预检

candidate2 manifest 绑定并哈希：`Rcode/multiFSYNC/R` 下全部 26 个 R 文件、DESCRIPTION、NAMESPACE、candidate2 runner/controller/evaluator/manifest builder/正式入口、五个实际 bayesSYNC reference R 文件、v0F 校准和 v0G/v0H/v0I/v0J 评价源码、两份协议。

另记录并验证：完整已安装 multiFSYNC namespace 函数 inventory/body/formals 哈希、公开主拟合/预得分函数哈希、R 版本、平台、BLAS、LAPACK 和关键依赖版本。formal 模式只接受状态为 frozen 且 `formal_execution_authorized=TRUE` 的另建绑定；源码、namespace、环境、seed/manifest 任一不一致立即拒绝。

## 12. 执行链、测试与产物

正式入口覆盖但本轮不运行：生成/密封 → B0/A/R/pooled → terminal → B0/A/R/pooled 真值盲选择 → offline racing → 03/08 continuation/horizon freeze → 授权解封 → multiFSYNC/pooled 评价 → cluster bootstrap/成功门汇总。

微型 smoke 必须在独立目录标记 `formal_experiment=FALSE`，执行同一顺序的缩小链，并明确不写 formal manifest 启动字段。回归测试至少覆盖公开/私有预得分等价、warning/error terminal、精确 join、formal fake/hash/seed 拒绝、四项 alignment、稳定性反例、truth lock、resume 幂等、audit error 和 pooled 全三方向/N/A/PPI 规则。

机器清单固定：10 数据、210 主 fit、72 audit records；seed collision audit 区分有意核心配对/continuation 原 seed 复用与非法碰撞。候选产物只标记 candidate，不可自行改为 frozen。

## 13. 时间预算与授权门

保守串行预留：200 条 multiFSYNC 主拟合约 61–67 小时（严格早停可降低实际值）；03/08 审计约 16.5–19.8 小时；10 条 pooled 约 10–20 小时；评价/哈希/I/O 约 3–7 小时，总计约 **91–114 小时**。8 路并发建议预留 **14–23 小时**墙钟，16 路约 **8–15 小时**，但内存和 I/O 可能限制缩放。

本文件完成后停在“candidate2 等待用户审阅”。只有用户先明确批准冻结，才可另建 frozen bundle；只有随后再次明确“开始正式实验”，才可生成十个正式数据或启动 210/72 条正式工作。

