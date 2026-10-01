# Stage 7E：Stage 7C 结构失败的定向机制诊断

日期：2026-09-28。范围：既有 Stage 7C 三个 C1 数据种子。性质：真值辅助初始化的事后机制诊断；不是可部署拟合方法、正式候选、独立 Monte Carlo 重复或论文性能分母。Stage 7C 的 108 条拟合、9 个真值盲 winner 与全部机器表均未回写。本轮未改模型、先验、CAVI、ELBO 或生成器。

## 1. 结论与证据边界

**在这三个固定数据集上，现有目标和原版 CAVI 可以到达结构正确、普通 ELBO 更高的局部端点；seed01/02 的 Stage 7C 失败主要定位为登记随机路径未产生这类候选。** 两个失败种子的真值投影起点与仅扰动三个活跃载荷方向的近真值起点，均保持完整共享—特异 core，且普通 ELBO 分别高于各自真值盲 winner 约 6000 和 8546。成功对照 seed03 的两条诊断端点也保持完整 core，比既有正确 winner 高约 201。六条均通过原版 Stage 7C 冻结结构评价口径、目标资格和严格实用收敛检查。

这证明的是**存在一个更高的正确局部盆地，且所测的狭窄载荷扰动能够留在其中**。不能由此声称正确解是全局最优、一般有限样本都可达、随机初始化只要增加次数就会找到，或某一条 PPI、尺度、更新公式已经被证实为根因。`near` 只在三个真实载荷的活跃支持内把方向余弦设为 0.995；均值、β、时间函数、得分、稀疏支持和噪声初值仍由真值辅助。因此它没有测量真实随机起点到正确盆地的体积。

## 2. 原始拟合与 FPCA 后处理

在三个数据种子的 C0/C1 ELBO winner，以及 seed01 C1 start11、seed02 C1 start08 两个 B2 方向活跃但 core 不完整的端点上，逐研究比较原始后验均值重构与 FPCA 后处理重构，共 8 个 fit、16 个研究响应。均值与 β 曲线最大绝对差均为 0；因子重构 RMS 差相对于原始因子 RMS 的最大值为 `5.389061e-05`，完整信号对应最大值为 `3.178386e-05`。逐载荷 PPI 最大差为 0，非零载荷方向的绝对余弦为 1。FPCA 规范化后的非零载荷范数/原始载荷范数在约 0.391–1.076 之间，所以载荷尺度指标必须注明使用原始还是规范化对象。

seed01/02 C1 winner 的研究 2 特异载荷在原始与后处理结果中，`PPI >= 0.5` 支持数都是 0；seed03 正确 winner 为 48。两个失败池各自存在 B2 局部活跃端点，支持数分别为 42 与 46，但另一结构块不完整。**B2 缺失已发生在拟合端点，不是后处理造成的。**直接机器依据：`results/RAW_POST_FACTOR.csv`、`results/RAW_POST_RECONSTRUCTION.csv`。

## 3. 六条定向拟合

每个数据种子有 `exact` 与 `near` 两条诊断 fit，共六条。初始时间函数和均值/β 投影到该数据绑定的 `K=9` O'Sullivan 基；真得分与载荷进入合法 continuation state，活跃/非活跃载荷 PPI 初值为 0.99/0.01，误差精度采用真值倒数，其余超参数继承 1-sweep 模板。每个种子的两条起点共用模板状态及拟合种子。六条均在 Stage 7C 同一绑定模型包上运行 400 次普通 `T=1` sweep，无退火；Stage 7C 正式随机起点为 99 次退火加 400 次普通 sweep。因此表中的终点普通 ELBO 可按同数据、同模型比较，而起点机制和计算路径不能作公平方法性能比较。

定义 `ΔELBO = 诊断端点普通 ELBO − 同数据 Stage 7C C1 真值盲 winner 普通 ELBO`。所有诊断 ELBO 均为目标资格有效值；无自适应 jitter 失效、无拟合 warning、无负的普通 ELBO 增量，400 次普通 sweep 均完成。

| 数据 | Stage 7C C1 winner：core / B2 | exact ΔELBO | near ΔELBO | 两条诊断：core / 精确因子结构 | B2 真支持找回数 | 最小载荷余弦 |
|---|---|---:|---:|---|---:|---:|
| seed01 | 否 / 关闭 | +6000.145 | +6000.061 | 均为是 / 是 | 42/50 | 0.99934 |
| seed02 | 否 / 关闭 | +8546.482 | +8546.404 | 均为是 / 是 | 46/50 | 0.99970 |
| seed03 | 是 / 保留 | +200.908 | +200.844 | 均为是 / 是 | 48/50 | 0.99960 |

三个 exact/near 配对的终点 ELBO 绝对差仅为 0.084、0.077、0.064；配对的 B2 真支持数相同。六条的共享、研究 1 特异、研究 2 特异三个方向均达到绝对载荷余弦 `>0.9993`；冻结评价器的全局正确计数 `3/3`、研究角色正确计数 `4/4`，缺失/错置/重复/额外因子计数为 0；真特征成分匹配 `8/8`。这里的“精确阈值因子结构”指因子数与角色位置，不代表 50 个变量的支持完全找回。直接机器依据：`results/oracle_summary_v1/evaluation_v2/ORACLE_ENDPOINT_COMPARISON.csv`、`ORACLE_FACTOR_DIAGNOSTICS.csv`、`ORACLE_FROZEN_EVALUATOR_SCIENCE.csv`；汇总在 `results/derived_v1/STAGE7C_ORACLE_CONTRAST.csv` 与 `ORACLE_PAIR_STABILITY.csv`。

上述 ΔELBO 是所有变分参数共同到达的端点目标差，不能拆解为 B2 单块的独立收益，也不能直接估计先验对 B2 的净惩罚。它足以反驳“现有错误 winner 已在正确结构之上”这一有限候选比较命题。

六条均标记 `strict_practical_converged=TRUE`；三个 Stage 7C C1 winner 也都标记为 TRUE，其中两个 core 错误。另按 Stage 7C runner 的完整定义 `practical_converged && converged && convergence_status=="converged_practical"` 复核六条，结果全部 TRUE，见 `results/STRICT_STATUS.csv`。因此该停止标签不是结构正确的充分条件；本实验也不能把它当作必要条件。六条在 sweep 380→400 的 ELBO 增益仅约 0.108–0.217；完整路径见 `results/derived_v1/ORACLE_ELBO_PATH_SUMMARY.csv`。在**固定 C1 维度配置**内，原真值盲最大合格 ELBO 规则没有机会选到 seed01/02 的正确结构候选；这些诊断端点若被允许进入池，ELBO 会高于原 winner，但它们实际永不进入正式候选池。此结论与跨 L/M 的维度选择无关。

## 4. 早期轨迹与尺度

诊断路径从第一次普通 sweep 开始记录 `1:20,25,30,40,50,75,100,150,200,300,380,400`。以 exact 路径的 B2 为例，第一次更新后平均载荷 PPI 为 seed01/02/03 的 0.0996/0.0980/0.1843，第 5 次为 0.0866/0.0933/0.0974，第 400 次为 0.0869/0.0936/0.0982。平均 PPI 约 0.1 与“500 个变量中约 50 个活跃”的稀疏结构相符，**不能单凭 0.1 判为关闭**；终点支持、载荷范数和方向也保持活跃。`near` 的轨迹与同种子 exact 很接近。直接轨迹为 `results/oracle_summary_v1/evaluation_v2/ORACLE_DRIVER_TRACE.csv`，精简快照为 `results/derived_v1/ORACLE_B2_TRACE_SNAPSHOT.csv`。

用连续网格真实时间函数 L2、真实得分 RMS 和真实载荷范数作分母，三个 exact 终点的 B2 原始载荷范数比为 1.085/1.059/1.089，时间函数 L2 比为 0.810/0.929/0.932，得分 RMS 比为 0.989/0.978/1.048。三个块确有尺度调节，seed01 的时间函数偏小，但六条均维持正确方向；没有出现仅属于两个失败种子的整体尺度坍缩。轨迹里的 `specific_time_function.posterior_mean_norm` 是**基系数欧氏范数**，不是连续时间函数 L2，因此不能用其大幅变化直接诊断函数能量爆炸。终点真函数尺度复核见 `results/ORACLE_B2_ENDPOINT_SCALES.csv`。

这些是**真值辅助路径**的轨迹。Stage 7C 随机路径只保存总迭代 1、380、400；总迭代 1 仍处于退火温度 1.9，不能与这里的第一次普通更新作同条件比较。随机路径在第 2–379 次究竟先丢方向、先降低 PPI，还是发生尺度反馈，仍未定位。

## 5. 竞争解释与下一分支

| 解释 | 本轮支持 | 反证或剩余限制 | 判定 |
|---|---|---|---|
| 有限样本 B2 真值共同过弱 | seed01 B2 观测 SNR 比 seed03 低；样本时点也略少 | seed02 B2 SNR 更高，Gram 与 K9 投影更有利；两失败种子均有更高 ELBO 的正确端点 | 不能解释两个失败池的共同现象 |
| 当前 ELBO/先验必然偏好缺 B2 | 原随机 winner 均关闭 B2 | 正确诊断端点在 seed01/02 的普通 ELBO 高约 6000/8546 | 对这三个已测池，不成立；不能推出全局排序 |
| 当前随机路径的联合盆地可达性不足 | 两失败池 12/12 无完整 core；六条真值辅助端点均完整，且目标更高 | 只测试狭窄真值邻域，随机路径的偏离时刻未知 | **当前主解释** |
| 早期稀疏门或尺度补偿触发失败 | 原随机池在后期有 B2 关闭；CAVI 无持续硬尺度规范 | 真值辅助路径不关闭；原随机路径缺中段记录，不能判断因果先后 | 待定，进入下一项定位 |
| FPCA 后处理或明显实现漂移 | 原始与规范化载荷尺度确有差异 | 原始端点已缺 B2；同一绑定包和冻结评价器可产出六个正确端点，输入/源码哈希及 K9 设计核验通过 | 不是当前首要解释；细微实现问题仍可查 |

**下一分支选“原随机路径的最小密集诊断重放”，而不是修改先验或 CAVI。**预先固定 Stage 7C C1 真值盲 winner 三条：seed01 start06、seed02 start11、seed03 start06；前两条为失败例，第三条为成功对照。使用原观测包、原种子、`random + pre_score=1 + G12 + anneal(1,1.9,100)`、相同 99+400 sweep 与全部超参数，只增加 `trace_sweeps = sort(unique(c(1:20, seq(25,90,5), 91:110, seq(115,200,5), 250,300,350,380,400,450,499)))`，覆盖退火前段、退火→普通转折和后段。真值不参与重放输入、起点或停止；完成后才用于解释方向。先用 1–2 次 sweep 微型 smoke 验证记录与包绑定，再启动最多三条；不自动追加起点。

重放的首要门槛是**端点复现**：普通 ELBO 差绝对值须不超过 `1e-6`，载荷 PPI/支持和原始载荷也须与原 Stage 7C 对应 fit 在数值容差内一致。若不一致，停止机制解释并审计诊断记录是否改变执行路径。若一致，输出逐总迭代的块级 PPI、载荷后验范数、得分范数、时间函数系数范数、温度与 driver/precision 比，定位 B2 后验质量第一次明显损失的区间，并比较 B1 与成功对照。现有 trace 不能恢复逐迭代载荷方向或真正的函数 L2；若 PPI 始终活跃而方向仍错误，或系数范数无法区分尺度反馈，应先设计仅供诊断的状态快照，再决定唯一一个真值盲干预。根据既有三条 Stage 7C winner 的实测拟合耗时约 4309–4526 秒/条，三条并行重放预计约 75–100 分钟墙钟时间，另需约 20–30 分钟准备、微型 smoke 与复核。该下一轮长拟合尚未启动。

该分支若获后续拟合授权，只新增 Stage 7F 独立 runner、协议和结果目录；不修改 `R/` 中的模型实现、任何更新公式或既有 Stage 7C/7E 文件。最小测试为 runner 解析、源/输入哈希、微型 smoke 与上述端点复现门槛；回滚边界就是新诊断目录，不借由 `reset`、`clean` 或覆盖既有结果回退。若需进一步记录逐迭代载荷向量或真正函数 L2，再单独列出诊断代码改动、公式与测试后实施。

## 6. 版本与复核

- 当前嵌套 Git 仓库：`dev/v0lv-init-geometry-20260818`，HEAD `5ec35267def4ccbc1a71c6811a9b06acccb19e9c`；保留全部既有 dirty/untracked 文件。Stage 7E 仅在自己的新实验目录写脚本、协议、报告和机器表。
- 本轮六端点解封评价的用户授权记录为 `STAGE7E_EVALUATION_AUTHORIZATION_20260928.txt`；真值从未进入 Stage 7C 的冻结候选选择。
- Stage 7C 冻结 marker 仍声明 `registered_fits=108`、`selected_winners=9`、`truth_used_for_selection=FALSE`；ECS 上 winner 文件 SHA256 `1e1e7791b3e975287ab195c87a801ab0fa24f97d26558a85fccc3b11af70a877`，与先前科学评价归档记录一致。
- Stage 7C 绑定模型源码归档 SHA256：`93b5fb86fc9cb687154502d00a0920979f51aa7f593966eabb462b2d73fcf924`。Stage 7E 六条拟合 runner SHA256：`35e0dab45bd4a1a68e05909c9ebdecd0eff8ae3775b03a0db1c664adcf4ff453`。拟合输入观测与已授权真值的逐数据哈希记录在 `results/oracle_summary_v1/INPUT_BINDING.csv`。
- 唯一执行的六端点评价脚本为 `evaluate_stage7e_oracle_20260928_v3.R`，本地与 ECS SHA256 均为 `72ece3a5db122567a04c090e44ae34703a5726a1f93ffa558899c9fab106ca2e`；冻结 Stage 7C evaluator 覆盖函数哈希也在 `evaluation_v2/EVALUATION_COMPLETE.txt`。`results/EVALUATOR_BINDING_V2.log` 记录它对既有 seed03 winner 的结构计数逐项复现。旧 runner 内的 `evaluate` 分支和独立评价脚本 v2 未用于本次六端点评价。
- 六条拟合主管状态 `6/6 TRUE`，六个 `FIT_COMPLETE.txt` 齐全、无 `FIT_ERROR.txt` 或 warning。只读尺度脚本 SHA256 为 `0b4704a52d9d6d563061e0266e15af0dfaa4e8461623fa874d06c9a5ecf9ed3a`；完整停止语义复核脚本 SHA256 为 `80385894791f957d253f7261f797f568004560cd6866a92db8d9150c6ad1cd12`；本地汇总脚本 SHA256 为 `46ae76098beec0e9fac117d370fc3da7847d81b39d571a3255b8ff9ca75c41ed`。
- 最终 ECS 汇总归档 `STAGE7E_SUMMARIES_20260928_v3.tar.gz` 的 ECS/本地 SHA256 均为 `60eead3a981feff0ca29f1ac2593ee4074ba65615b7f742612ca099c50f5cecc`，63 个条目，无 fit RDS 或真值 RDS。原始 fit RDS 和 sealed truth 保留在 ECS 原位置；先前 v1/v2 归档保留，v3 增加了完整停止语义核验。
- 六个原始 fit RDS 与六个初始化 state RDS 另以 `STAGE7E_FITS_AND_STATES_20260928_v1.tar.gz` 回传本地；ECS/本地 SHA256 均为 `95b20a059f10345971e18526db72df7c2857f2789ef47051888849e9384798e7`。此归档不含 Stage 7C 真值包，但 state RDS 包含真值辅助初值，只能用于机制诊断，不能进入真值盲选；原始文件继续保留在 ECS。

理论边界：这些结果针对生成器命名的三个 rank-one 载荷方向和固定 `L_f=1,L_s=(1,1),M=2,K=9`。它们不补足一般共享—特异因子的唯一识别证明；DGP 的载荷列平方范数相等，也不能据此声称满足文档中的严格递减强识别条件。未知维度的 L/M 选择仍是独立问题。
