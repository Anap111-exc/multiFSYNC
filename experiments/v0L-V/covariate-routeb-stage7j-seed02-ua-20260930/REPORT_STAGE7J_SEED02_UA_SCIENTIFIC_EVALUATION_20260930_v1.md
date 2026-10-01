# Stage 7J：seed02/C1/start11 幅值恢复复制的科学评价

日期：2026-09-30。范围：Stage 7C 已冻结的 `covrb7c_02__C1__start11` 从 sweep 50 的一次 U/A 配对续跑；这是同一个数据集上的诊断分支，不是新增独立数据重复或正式候选。真值仅在盲层端点冻结、用户授权后用于评价。模型、先验、CAVI、普通 ELBO、更新顺序、温度表、PPI 规则及固定维度均未修改。完整原始 fit 和状态保留在 ECS `/root/multifsync-results/`；本机只保存小型机器表。

## 1. 协议与可比性门槛

预登记协议见 [PROTOCOL_STAGE7J_SEED02_UA_20260930_v1.md](PROTOCOL_STAGE7J_SEED02_UA_20260930_v1.md)。U 从原状态原样续跑；A 保留原 B2 载荷、得分和时间函数的均值方向，只按 Stage 7I 已固定的规则恢复幅值。A 的产品均值 RMS 与预构造的真值盲 S 状态同为 `0.2907836345`，有效载荷范数同为 `7.87522814`；S **只用于幅值标定，没有拟合**。保留原 PPI、slab 方差、得分与函数协方差、其他块状态；均值变化使载荷二阶矩和后续精度也变化，因此 A 不是单独固定一切中介量的“纯尺度”试验。

U、A 各完成 sweep 51–499：49 次剩余退火加 400 次普通更新，共 898 次后续更新。U 对冻结 Stage 7G 的 449 个状态快照作 19,307 项字段比较，最大状态差和普通 ELBO 差均为 0，PPI 零模式一致；因此续跑接入门槛通过。盲层端点先冻结；后续真值评价的 U 与原 Stage 7C start11 科学表 231 个字段均一致，最大数值绝对差 `4.075×10⁻¹⁰`。U 与 A 的目标资格、停止状态又以原口径单独复核，均 `objective_eligible=TRUE`、`strict_practical_converged=TRUE`、400 次普通更新完成、无普通 ELBO 下降或 T=1 solver jitter。严格实用收敛只描述停止语义，不证明结构正确。

直接证据：[盲层端点冻结表](stage7j_small_readout_20260930_v1/stage7j_seed02_blind_freeze_20260930_v1/BLIND_ENDPOINT_SUMMARY.csv)、[停止语义复核](stage7j_small_readout_20260930_v1/stage7j_seed02_status_semantics_20260930_v1/OBJECTIVE_AND_STOP_STATUS.csv)、[真值评价](stage7j_small_readout_20260930_v1/stage7j_seed02_scientific_evaluation_20260930_v1/STAGE7J_ENDPOINT_SCIENCE.csv)、[U 科学评价一致性](stage7j_small_readout_20260930_v1/stage7j_seed02_scientific_evaluation_20260930_v1/U_BASELINE_SCIENCE_PARITY.csv)。停止语义表 SHA256 为 `afb76a205fc0bb8e5d2e128027d8199dd61c67283ca5063450514f11d438c7b7`；科学表 SHA256 为 `abb8f988377e333053b0cdfc6cf38c52c0103abb306c149435ce63d22bbd994a`。

科学表中的 A `strict_practical_converged` 留空，因为通用真值评价器没有为新诊断分支重算该元数据；上述单独的无真值停止语义复核读取了冻结 fit 的完整状态，给出 A 为 TRUE。以停止语义复核表为准，保留原科学表不回写。

## 2. 端点结果

| 固定 seed02/C1/start11 指标 | U 原状态续跑 | A 幅值恢复 |
|---|---:|---:|
| 普通 ELBO | −173563.471994 | −173564.467705 |
| objective eligible / strict practical converged | 是 / 是 | 是 / 是 |
| 完整 core / 精确阈值结构 | 否 / 否 | 否 / 否 |
| shared / B1 / B2 保留数 | 1 / 1 / 0 | 1 / 1 / 0 |
| 真函数分量覆盖 | 6/8 | 6/8 |
| 观测点完整信号 NRMSE | 0.363207327 | 0.363207381 |
| 稠密网格完整信号 NRMSE | 0.395326659 | 0.395328981 |
| 完整角色贡献平均 NRMSE，缺失计零 | 0.367319456 | 0.367329048 |
| β 相对 ISE | 0.078588866 | 0.078589094 |

A−U 的普通 ELBO 为 `−0.995711`，观测点完整信号 NRMSE 差仅 `+5.37×10⁻⁸`。两支都正确放置 shared 与 B1 的研究角色，研究 2 的真实 B2 角色均缺失；shared/B1 载荷方向几乎不变，未见其他块的实质端点损伤。`R_A/R_B1/R_B2`：U 为 `0.999783814/0.999382502/0`，A 为 `0.999783863/0.999382759/0`。这些结论来自端点真值评价，不是根据 ELBO、PPI 或停止标签推测。角色明细见 [STAGE7J_ROLE_DAMAGE.csv](stage7j_small_readout_20260930_v1/stage7j_seed02_scientific_evaluation_20260930_v1/STAGE7J_ROLE_DAMAGE.csv)，成对差见 [STAGE7J_PAIRED_DIFFERENCES.csv](stage7j_small_readout_20260930_v1/stage7j_seed02_scientific_evaluation_20260930_v1/STAGE7J_PAIRED_DIFFERENCES.csv)。

## 3. 已保存轨迹中的第一个可见分叉

初始 A 的 B2 均值产品 RMS 为 `0.290784`，完整后验产品二阶 RMS 为 `0.723250`，但 PPI≥0.5 的载荷个数仍为 0。这是注入后的状态，不是一个已完成的坐标更新。第一次续跑后：

| 总 sweep | A 有效载荷范数 | A 均值产品 RMS | A 完整后验二阶 RMS | A 的 B2 PPI≥0.5 数 |
|---:|---:|---:|---:|---:|
| 50，注入后 | 7.8752 | 0.290784 | 0.723250 | 0 |
| 51 | 19.3494 | 0.019963 | 0.026478 | 125 |
| 52 | 7.3820 | 0.006712 | 0.014183 | 35 |
| 53 | 4.5245 | 0.003734 | 0.011814 | 11 |
| 55 | 2.4335 | 0.001643 | 0.009491 | 1 |
| 56 | 1.8713 | 0.001159 | 0.008757 | 0 |
| 60 | 0.02535 | 2.729×10⁻⁷ | 0.005337 | 0 |
| 65 | 0 | 0 | 0.003701 | 0 |
| 499 | 0 | 0 | 5.652×10⁻⁷ | 0 |

因此，sweep 51 的初始产品缩小**不能简单归因于载荷或 PPI 立即关闭**：有效载荷范数反而增加，PPI 支持暂时升到 125；同时得分—时间函数均值乘积范数从 `12.1625` 降到 `0.3398`。之后支持数与产品共同下降；到 sweep 65 均值为零，但完整二阶矩仍非零，不能说整个 B2 后验在该时刻已经消失。已存快照只给出每个完整 sweep 的状态，不能单靠这张表确证某一个内部更新步骤是因果起点。原始轨迹见 [STAGE7J_B2_TRAJECTORY_SELECTED.csv](stage7j_small_readout_20260930_v1/stage7j_seed02_blind_trajectory_readout_20260930_v1/STAGE7J_B2_TRAJECTORY_SELECTED.csv)。

相同固定 A 规则在 Stage 7I seed01/C1/start06 中，sweep 51 也从 `0.136908` 降至 `0.004503`，但 sweep 52 回升至 `0.012293`、sweep 53 至 `0.070608`、sweep 55 至 `0.197855`，终点为 `0.257763`，B2 保留 42 个支持并恢复完整 core，A 相对 U 的 ELBO 高 `5894.132`，观测点完整信号 NRMSE 由 `0.286750` 降至 `0.141321`。见 [Stage 7I 科学报告](../covariate-routeb-stage7i-uas-20260929/REPORT_STAGE7I_UAS_SCIENTIFIC_EVALUATION_20260929_v1.md)。两组都是**已知失败病例上的一次诊断**；这一对比说明首次缩小不是必败判据，幅值注入也不是这两个病例都有效的固定修复。

## 4. 方向与竞争解释

真值只在冻结后用于解释注入方向。初始 A 的有效 B2 载荷与真实 B2 的绝对余弦：seed01 为 `0.332015`，seed02 为 `0.044095`；seed02 A 与真实 shared 的绝对余弦为 `0.313039`。seed02 A 在 sweep 51 与真实 B2 的余弦进一步降至 `0.012361`；到 sweep 60 虽升至 `0.243444`，有效载荷范数已仅 `0.02535`，随后归零。真值盲 S 构造方向与真实 B2 的余弦为 `0.998401`，但 S 未续跑，**不能推断 S 会救回 seed02**。初始 U 的载荷极小，其方向余弦只能描述数值方向，不能代表已有可用信号。直接表：[LOADING_DIRECTION_TRUTH_DIAGNOSTIC.csv](stage7j_small_readout_20260930_v1/stage7j_seed02_initial_loading_truth_diagnostic_20260930_v1/LOADING_DIRECTION_TRUTH_DIAGNOSTIC.csv)。

seed02 的 A 幅值恢复失败，同时其初始方向比 seed01 更偏离真 B2，支持“方向／时间组织与反馈共同决定能否反弹”这个竞争解释，不能把方向差判成已经隔离的唯一原因。A 增幅改变二阶矩与随后的精度、PPI；残差首变量方向可见并不保证现有得分—函数产品沿真实 B2 形成有符号、可持续的坐标驱动。固定残差方向上的 signed-cancellation 比率只刻画该方向，不能等同于完整更新的抵消。载荷 trace 的 11 维 `data_driver_norm` 与 2 维 `linear_driver_norm` 不在同一空间，不能相减当成投影损失。现存逐 sweep 快照也不足以独立验证 `Q^{-1}h` 的各个中介环节。

有限样本真 B2 过弱不是这两个失败池的共同充分解释：此前 Stage 7D 发现 seed02 的 B2 观测信号 RMS 高于 seed01/03，Gram 与 K9 投影条件也更有利；Stage 7E 的真值辅助 seed02 诊断端点保持正确结构，普通 ELBO 比原真值盲 winner 高约 `8546.48`。这些端点**不属于正式候选池**，不能拿来计算可部署方法的成功率，也不证明随机路径有足够概率进入正确盆地。参见 [Stage 7D 报告](../covariate-routeb-stage7d-readonly-attribution-20260928/REPORT_STAGE7D_READONLY_ATTRIBUTION_20260928_v1.md) 与 [Stage 7E 报告](../covariate-routeb-stage7e-mechanism-20260928/REPORT_STAGE7E_MECHANISM_ATTRIBUTION_20260928_v1.md)。

## 5. 分支判断

预登记的 seed02 U→A 复制**失败**。目前没有根据把一次性幅值恢复直接并入正式方法，也没有根据立即修改先验、PPI、shared 更新顺序，或继续追加谱重启、幅值网格及密集重放。原模型和 CAVI 可以保留为研究基线；论文必须如实限定固定维度下的随机路径结构稳定性。若后续要争取更稳的结构恢复，应先预先规定一个真值盲干预与独立数据验证，和当前方法在相同预算下配对比较。若不开展该项，则报告已有失败及其限制，不能只挑 seed01 的成功干预或把 seed02 失败简单归为困难数据种子。固定 C1 起点可达性、跨 L/M 维度选择、一般识别性和 strict convergence 是不同问题。

评价代码：[evaluate_stage7j_seed02_endpoints_20260930_v2.R](evaluate_stage7j_seed02_endpoints_20260930_v2.R)（SHA256 `d249893de9be558a0cebc1ba204fc8c575b310e80a600eb8c0e41e53b2e48e85`），[analyze_stage7j_seed02_trajectory_20260930_v1.R](analyze_stage7j_seed02_trajectory_20260930_v1.R)（SHA256 `2f5d27144cf45c40ed19a8dec3ac8e0032b190a4e015fec5704cadff1ba48edd`），[evaluate_stage7j_initial_loading_alignment_20260930_v1.R](evaluate_stage7j_initial_loading_alignment_20260930_v1.R)（SHA256 `9e1b5fc621c711748f8816051dbe8fe6bea3a87ba1536037ab154460c09db433`），[audit_stage7j_seed02_status_20260930_v1.R](audit_stage7j_seed02_status_20260930_v1.R)（SHA256 `62e0f096b433844ba05b3fafcf919e0795d720b45337537140a802108e1c0bd6`）。较早的 `evaluate_stage7j_seed02_endpoints_20260930_v1.R` 因真值哈希字面量末位错误而在打开真值前的 preflight 停止；v2 修正并完成，未修改拟合端点。
