请独立审计 multiFSYNC Stage 7K 的 seed02/C1/start11 诊断，并提出下一项最小研究任务。先只读；不要新增 fit、修改模型/先验/CAVI/ELBO、重选端点、提交 Git 或停机。ECS 是用户为查看结果手动重新启动的。

项目根目录：`D:\文档\Factor model\multiFSYNC`。本地协议、代码、报告及约 0.5 MB 小表在 `Rcode\multiFSYNC\experiments\v0L-V\covariate-routeb-stage7k-seed02-s-20260930\`，请先读 `PROTOCOL_STAGE7K_SEED02_S_20260930_v1.md`、`REPORT_STAGE7K_SEED02_S_SCIENTIFIC_EVALUATION_20260930_v1.md`，再直接核对 `stage7k_small_readout_20260930_v1/` 的 CSV、哈希和 marker；报告不是高于机器表的证据。

ECS 可用 Workbench 连接：实例 `i-bp1egmsairexm9w1rnks`、地域 `cn-hangzhou`，当前公网 IP `47.98.87.163`（IP 可能变化，请优先按实例 ID）。执行代码 `/root/stage7k-exec-20260930-v1/`；补做的只读轨迹代码 `/root/stage7k-readout-20260930-v1/analyze_stage7k_s_trajectory_20260930_v1.R`。原始 S fit/快照 `/root/multifsync-results/stage7k_seed02_s_20260930_v1/S_full_v1/`；盲层冻结 `/root/multifsync-results/stage7k_seed02_blind_freeze_20260930_v1/`；真值评价 `/root/multifsync-results/stage7k_seed02_scientific_evaluation_20260930_v1/`；只读轨迹 `/root/multifsync-results/stage7k_seed02_trajectory_readout_20260930_v1/`。U/A 沿用 Stage 7J 已冻结端点，没有重跑。

设计：S 初态是 Stage 7J 解封真值前由观测残差构造并冻结的 sweep-50 状态，SHA256 `72a19fd46ed9f096966eca6c5b6f9bd673a6010e51a356ddcaa7e5156ad6d292`；**决定追加 S 是看过 seed02 真值之后**，故这是探索性病例诊断，不是前瞻性真值盲成功率。只续跑 51–499，共 449 sweep，原模型/先验/更新顺序/预算；先冻结 U/A/S 端点，再做同一 Stage 7C 口径科学评价。A 与 S 的初始 B2 均值产品 RMS 均 `0.290784`，但完整二阶矩、方向及时间/个体组织不同。

已核实结果：监督器 14:07:47 正常完成，S 耗时 6951.847 秒；U 对 Stage 7C 科学基线 231 项、U/A 对 Stage 7J 151 项比较均通过。U/A 的 B2 R=0、块数 1,1,0；S 的 B2 R=`0.999608`、块数 1,1,1，核心/精确块结构/8 个真特征覆盖均为 TRUE，shared/B1 R≈`0.999784/0.999426`。普通 ELBO U/A/S=`−173563.472/−173564.468/−165125.693`，三者 objective eligible 与已保存 strict practical 标签均为 TRUE。观测点完整信号 NRMSE U→S `0.363207→0.130466`，密集网格 `0.395327→0.204122`；S 不是零误差或完美逐变量支持。S 的 400 个普通 ELBO 点单调增加。

轨迹更关键：A/S 在 sweep 50 的 B2 均值产品 RMS 同为 `0.290784`；sweep 51 A/S 分别降至 `0.019963/0.004240`，然后 sweep 53 为 `0.003734/0.118451`，sweep 55 为 `0.001643/0.287254`，终点为 `0/0.363710`。A 的 B2 PPI≥0.5 数从 sweep 51 的 125 降至 60 的 0；S 从 51 的 19 增至 53 的 48，终点 46。即两支首次都缩小，随后是否再增长发生分叉；不能把高初始幅值或第一步 PPI 峰值等同于持续恢复。

请完成：
1. 先独立核对冻结哈希、端点科学表、轨迹与上述数字；指出任何口径或实现问题。
2. 判断本实验能支持哪些机制说法，尤其“seed02 的 S 整套状态比固定 A 规则更有效”与“谱载荷方向单独起效”之间的证据界限。检查初始二阶矩、时间函数、得分、载荷、PPI、后续精度与 shared/B1 变化；若需要，可在 ECS 对已存 fit/快照做只读分析。
3. 对比 seed01 Stage 7I（A/S 都成功）和 seed02 Stage 7J/7K（A 失败、S 成功），说明哪些结论可以跨两病例提出，哪些不能。不要把病例选择、strict convergence、固定 C1 起点救援与跨 L/M 维度选择混同。
4. 提出一项**前瞻、真值盲、最小**的下一实验：固定 S 构造/触发规则、输入、对照、计算预算、机器输出、科学成功标准及停止规则。若证据不足以改算法，也请直说。不要根据这两个已知失败 seed 继续回调谱方向或扩大无方向起点搜索。
