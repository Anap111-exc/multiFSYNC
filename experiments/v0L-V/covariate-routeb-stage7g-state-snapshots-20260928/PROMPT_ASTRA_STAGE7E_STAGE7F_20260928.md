请作为 multiFSYNC 的独立统计方法审计者，复核以下两轮机制实验。项目根目录是 `D:\文档\Factor model\multiFSYNC`；先读 Stage 7E、7F 两份报告及其引用的机器表，必要时再追溯冻结代码。交接和下述文字只是导航；以冻结输入、哈希、runner 和机器结果为准。

背景：Stage 7C 固定 C1 维度的真值盲候选中，seed01/02 各 12 个起点均没有正确核心，主要缺研究 2 特异块 B2；seed03 有正确 winner。拟合选择没有使用真值或结构标签。

**Stage 7E（6 条真值辅助机制诊断 fit）：** 对同三份数据分别做投影真值 exact 与仅扰动三个活跃载荷方向的 near 起点。seed01/02 的两类诊断端点均保留完整正确核心，普通 ELBO 相对各自 Stage 7C C1 winner 高约 6000.145/8546.482；seed03 高约 200.908。六条都通过客观 ELBO 资格和严格实用收敛检查。原始载荷与 FPCA 后处理比较也确认 B2 缺失发生在拟合端点，不是后处理删除。结论仅限：已测数据及当前目标/CAVI 中存在更高 ELBO 的正确局部盆地，窄真值邻域可保持；不证明全局最优、随机可达或先验/某一步更新的独立因果作用。Oracle 路径绝不进入正式候选池。

**Stage 7F（三条冻结 C1 真值盲 winner 密集重放）：** seed01 start06、seed02 start11、seed03 start06 使用原观测、原起点种子、原 99 次退火+400 次 T=1、原模型和 G12，只增加轨迹记录。三条完整 400 点 ELBO、54 项原始参数、9 组旧 trace 与 Stage 7C 均为零差，零 fit warning。失败例 B2 的得分/有效载荷均值在早期已经很弱，成功对照也曾短暂变弱但在约 sweep25–60 恢复。到 sweep50 三条 B2 平均 PPI 仍几乎一样（约 0.220–0.223），但有效载荷范数分别为 0.000692、0.0000723、0.811。seed01/02 的 PPI 到约 sweep65 跌至 <0.01，约 sweep70（温度约 1.215）联合关闭并持续到终点；seed03 在 55–60 出现支持与载荷/得分同步回升。代码中得分更新依赖当前有效载荷、载荷更新依赖得分和时间函数，支持联合反馈后由稀疏门控锁定；“首 sweep 由 PPI 直接关闭”与轨迹不符。当前 trace 不能恢复逐 sweep 载荷方向或真实函数 L2，因而还无法判断起点方向从未正确、方向先丢失、或尺度先坍缩的具体顺序。B1 轨迹也不支持简单的 B1 单独抢占解释。

两份报告：
- `Rcode/multiFSYNC/experiments/v0L-V/covariate-routeb-stage7e-mechanism-20260928/REPORT_STAGE7E_MECHANISM_ATTRIBUTION_20260928_v1.md`
- `Rcode/multiFSYNC/experiments/v0L-V/covariate-routeb-stage7f-dense-replay-20260928/REPORT_STAGE7F_DENSE_REPLAY_ATTRIBUTION_20260928_v1.md`

请明确区分观测证据、机制推断和未证实假设；指出任何与机器表或代码不符之处；给出最有区分力的下一项最小诊断及条件化的单一真值盲算法干预、比较与停止标准。当前 Stage 7G 已获授权重放这三条路径，仅增加每 sweep 原始状态快照；必须先验证端点及 Stage 7F 轨迹一致，再用已授权解封的真值做方向评价。请勿用真值、结构标签、NRMSE、ISE 或 R/P/L 选择拟合；不要把 strict convergence 当作结构正确的充分/必要条件，也不要把固定维度起点选择混同跨 L/M 维度选择。
