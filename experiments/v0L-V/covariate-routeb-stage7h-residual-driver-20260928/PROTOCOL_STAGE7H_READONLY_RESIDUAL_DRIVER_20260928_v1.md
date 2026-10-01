# Stage 7H 只读残差与更新驱动归因协议

日期：2026-09-28。授权范围：按用户对 ASTRA 审计的指示，利用 Stage 7G 已冻结三条 C1 winner 的观测包、fit、全状态快照与 driver trace 做一次只读归因。**不调用 `bayesSYNC_multi()`，不运行新 fit，不改变模型或已有结果，不产生正式候选。** 原始真值只在真值盲表、脚本哈希和表哈希绑定之后用于独立评价。

## 问题与事先界定

Stage 7G 逐 sweep 真值对齐表表明：成功 seed03 的 B2 在 sweep 60 与真共享载荷方向余弦 `0.998726`，而与真 B2 余弦 `0.000659`；拟合共享载荷同时与真共享方向余弦 `0.998657`。seed01/02 的拟合共享载荷在 sweep 40 已较早对齐真共享方向（`0.99918/0.99987`）。待判别的是研究 2 的共享方向残差信号是否促成 B2 暂时激活，以及失败路径中真 B2 残差信息是否存在却未被当前更新利用。载荷余弦只说明空间方向，不等同于完整时间函数、个体得分或信号重构。

对研究 2、受试者 `i`，B2 留出残差定义为工作尺度的

`R_i^(-B2) = Y_i - C_i nu_mu - Σ_r Z_ir C_i nu_beta,r - a ⊙ (C_i nu_phi zeta_i)`。

该研究只有 B2 一个特异块；研究 1 的 B1 不在该残差中。所有拟合观测计算保持原观测时间、原 `C_i`、原似然精度与原变量顺序。另列不扣除均值、β、共享项的残差版本，用于定位扣除哪一项时可见方向改变。这是固定端点状态上的代数比较，不代表重新拟合的因果效应。

## 输入、时点与冻结顺序

- 输入绑定：Stage 7G 的 [FIT_MANIFEST.csv](../covariate-routeb-stage7g-state-snapshots-20260928/stage7g_20260928_v1/FIT_MANIFEST.csv)、[INPUT_BINDING.csv](../covariate-routeb-stage7g-state-snapshots-20260928/stage7g_20260928_v1/analysis_v1/INPUT_BINDING.csv)、三个 `RAW_STATE_SNAPSHOTS.rds`、三个已绑定的 Stage 7C observation bundle、Stage 7G 冻结源码和 Stage 7F/7G driver trace。读取前重核每个输入的 SHA256、winner ID、迭代索引 `0:499`；哈希不符即停。
- 主分析检查 `t=20,25,40,50,55,60`，三条路径各读取 `t-1` 与 `t` 共 18 个更新时点。关闭窗口仅 seed01/02 增加 `66→67`。先以 seed03 `49→50` 做一个重建点的微型 smoke，再决定是否跑完整固定范围。每次只加载一条路径，完成后释放对象。
- 真值盲层先输出 `RECONSTRUCTION_PARITY.csv`、`RESIDUAL_GEOMETRY_BLIND.csv`、`UPDATE_DRIVER_BLIND.csv`、`BLOCK_ENERGY_BLIND.csv` 与文件 SHA256。全部检查通过并锁定后，另读已授权真值，输出 `RESIDUAL_TRUTH_ALIGNMENT.csv`。不得根据真值改变时点、分组、谱秩或盲层公式。

## 重建公式与核验

为复现 sweep `t` 内的 B2 载荷更新，均值、β、共享/特异时间函数、共享/特异得分和共享载荷取 `t` 的快照；噪声精度与 B2 的 omega 对数矩取 `t-1` 快照，因为它们在 B2 载荷之后更新。退火系数 `c_t=1/T_t` 取 `t`。不可将 sweep 末所有字段直接代入 sweep 内 driver trace。先按源码逐项复现一个 `t=50` 的 `data_driver_norm`、`precision_trace`、`slab_mean_norm`、PPI 范围/均值及快照中的载荷和 PPI。

每坐标精度与驱动分解为

`h_j = tau_j Σ_i gbar_i^T R_ij^(-B2)`, `U = Σ_i E||g_i||²`, `v_j = {c_t(tau_j U + 1)}^{-1}`, `m_j=c_t v_j h_j`, `bbar_j=rho_j m_j`。

`gbar_i = C_i nu_psi xi_i`，`E||g_i||²` 同时包含函数/得分均值及二阶协方差。PPI logit 必须包含冻结 `update_loadings.R` 的温度项 `0.5(1-c_t)log(2π)`，并区分 prior log odds、`m_j²/v_j` 与 `log(v_j)`。记录每坐标 `h_j,v_j,m_j,rho_j`，求和与相消程度，记录当前函数空间捕获率及方差项对 `U` 的贡献。得分线性项和精度按同一相邻状态顺序重建；在任一重建未通过前，不解释其因果含义。

残差几何分别在原观测残差、及由原 `C_i` 投影到拟合基的残差上计算变量维度的前 5 个奇异方向/值和谱间隙。个体分组固定为受试者索引奇/偶；投影所用 QR 阈值固定 `1e-8`，记录每个设计矩阵的有效秩、时点数和条件数。若 `C_i` 已有满行秩，则观测点投影为恒等算子，这一“基内捕获率=1”没有诊断价值，必须明确标注。若某一方向谱间隙不充分，只报告前 5 维子空间，不解释单根向量。噪声白化若用于谱计算，真值方向必须经过相同变量权重变换后才对齐。`BLOCK_ENERGY_BLIND.csv` 计算 shared、B2 的均值重构 RMS、二者内积与完整因子产品二阶 RMS；逐 sweep 的旧 `conditional_second_signal` 不含时间函数协方差，不可充作完整二阶量。

## 输出与判定

| 表 | 最小字段 |
| --- | --- |
| `RECONSTRUCTION_PARITY.csv` | 输入哈希、fit ID、时点、观测设计及残差恒等式误差、driver/precision/slab/PPI 与保存 trace/快照的绝对和相对差 |
| `RESIDUAL_GEOMETRY_BLIND.csv` | 原观测及基内残差的前 5 个奇异值、谱间隙、奇/偶组前 1/2/5 维稳定性、基投影捕获率、设计时点数/秩/条件数 |
| `UPDATE_DRIVER_BLIND.csv` | 每坐标 `h_j,v_j,m_j,rho_j`、logit 分项、得分驱动/精度项、方差贡献、相消比和小项跳过数 |
| `BLOCK_ENERGY_BLIND.csv` | shared/B2 均值贡献、完整二阶产品 RMS、两个重构内积/余弦、残差扣除阶段的范数变化 |
| `RESIDUAL_TRUTH_ALIGNMENT.csv` | 盲表绑定后，前 5 维残差子空间对真共享/真 B2 的投影平方、两种方向在各扣除阶段的变化 |

“残差仍有可用 B2”的筛查证据预设为：真 B2 对前 5 维残差子空间的投影平方至少 `0.5`，并在两个相邻检查点及固定奇/偶组内稳定出现。该阈值不是显著性检验；低于阈值不能证明信息不存在。若残差有 B2 而基捕获/标准化驱动低，优先审计投影与更新几何；若驱动充足而方差精度或 PPI logit 使有效均值衰减，审计反馈；若真 B2 只在扣除某项后消失，审计该项重构竞争。成功例若 B2 的再增长驱动沿共享方向，需要把角色暂时重复与其后迁移当成独立机制分析。证据不分离时结论保持未决。

## 资源、失败与回滚边界

微型 smoke 封顶 `30 CPU 分钟`、`6 GB` 单进程内存；超出即保留日志并停止扩展。活跃量的重建核验容差事先固定为 `|差| ≤ 1e-7 + 1e-7×max(1,|参考值|)`；零值、极小驱动和阈值跳过逐项另核。哈希/设计不符、非有限值或复现失败时，不产出机制结论，也不启动拟合。新增文件均置于本 Stage 7H 目录，原 Stage 7C–7G 证据与项目 `R/` 均保持原状；放弃本诊断只需忽略本目录，无需 reset、clean 或覆盖任何用户文件。
