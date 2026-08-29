# G12 机制审计阶段 4F 完整报告

报告日期：2026-08-28
性质：对已完成 Stage 4E 的只读、开发性、已解封真值机制审计
正式 v0L-V 结果：否

## 1. 执行完整性

- 输入为 Stage 4E 已冻结的 24 个 fixed-400 端点：3 个数据集、每个 8 个起点；
- 3 个 winner 继续使用解封前按最大有效普通 T=1 ELBO 冻结的身份；
- 24/24 fit、24/24 checkpoint 对象和 24/24 真值评价完整；
- 1440 条 block-PPI checkpoint、1512 条 block-scale checkpoint、504 条全局 scale checkpoint、96 条 study-role 过程/贡献记录全部生成；
- 13/13 审计 QC 通过；运行耗时 32.21 秒，峰值内存约 244 MiB；
- 没有生成数据、重新拟合、continuation、重新选择 winner 或修改历史对象。

## 2. sparse_01：错误发生在很早的块竞争阶段

`g12ss4c_sparse_01` 的 8 个起点均未达到完整正确结构。checkpoint 表明，这不是 200–400 sweep 之间逐渐形成的错误：

- 7/8 起点的研究 1 特异块 B1 在第一个已记录普通 T=1 checkpoint（sweep 20）时，factor PPI 已降到约 `5.8e-8`–`6.7e-8`，有效载荷、feature、score 和贡献均为 0；
- 唯一未关闭 B1 的 seed 7 在 sweep 20 前关闭了共享块 A，factor PPI 约 `3.7e-8`；
- seed 2 和 seed 4 还同时关闭了 B2；
- winner seed 6 的 B1 factor PPI 从 sweep 20 的 `5.8887e-8` 到 sweep 400 的 `5.8250e-8`，没有恢复迹象；A 与 B2 则始终为 1。

初始化 checkpoint 0 时，三个块的载荷均非零，范数约为真值的 1.5–1.8 倍；到普通 sweep 20，失败块已精确归零。由于 Stage 4E 没有记录 99 个退火 sweep 内部 checkpoint，目前只能把塌缩定位为“初始化之后、退火或最初 20 个普通 sweep 内”，不能进一步声称发生在某一个具体更新。

因此，`sparse_01` 的主要问题是早期共享—特异块竞争/关闭，而不是 400 sweep 不足，也不是把最终 PPI 阈值从 0.5 调成 0.6–0.9能够解决的问题。继续 winner 到 800 没有针对性价值。

## 3. sparse_02：身份正确与科学尺度正确是两件事

Stage 4E 的 sparse_02 winner 为 seed 7，三个 factor PPI 在所有记录 checkpoint 均为 1，三个同角色载荷方向在 sweep 400 的余弦为：

- A：0.99954；
- B1：0.99926；
- B2：0.99916。

原始 CAVI 状态中的载荷范数比在 sweep 400 并不爆炸，分别为 1.457、1.100 和 1.217；但最终正交化/报告尺度下的载荷范数比为 7.009、3.669 和 0.932。与此同时，dense factor-process 已在持续放大：

| 量 | sweep 200 | sweep 400 | 变化 |
|---|---:|---:|---:|
| A feature grid RMS | 1.473 | 3.162 | 2.15 倍 |
| A factor-process grid RMS | 2.286 | 5.465 | 2.39 倍 |
| A complete-contribution grid RMS proxy | 1.101 | 2.389 | 2.17 倍 |
| B1 feature grid RMS | 2.022 | 2.654 | 1.31 倍 |
| B1 factor-process grid RMS | 2.976 | 3.644 | 1.22 倍 |
| B2 complete-contribution grid RMS proxy | 0.303 | 0.285 | 基本稳定 |

A 的过程和贡献在 340→400 仍单调增加；B1 在中途发生明显尺度重分配后仍未完全稳定。第 400 步的 short RSS 相对变化为 0.00362，long fitted NRMSE 为 0.00511，long RSS 相对变化为 0.03587，分别超过 0.001、0.003、0.006 门槛。

这解释了为何该端点能同时出现：

- 观测点 signal NRMSE 仅 0.211；
- 观测时点完整贡献 NRMSE 均值为 0.311；
- 但 dense-grid signal NRMSE 达 3.580，甚至远差于 mean-only dense NRMSE 0.702；
- covariance operator relative error 达 40.13。

结论不是“载荷方向错了”，而是稀疏时间观测下，正确身份的功能成分在观测点之外发生严重的尺度/形状失控。代码复核表明，最终 FPCA 后处理使用完整后验二阶矩得到 `factor_scale`，把载荷乘以该尺度、把得分除以同一尺度；因此它按设计保持完整贡献，而不是凭空制造大贡献。7.009/3.669 的 canonical loading norm ratio 主要是在报告参数中显式暴露了已经膨胀的因子过程积分方差。400→800 审计用于判断这种底层过程尺度是否仍在演化，而不是把问题简单归咎于一次符号或排列后处理。

## 4. 补齐的因子过程与完整贡献结果

三个 truth-free winner 的 loading-first matched 结果为：

| 数据 | 结构 | matched/total | factor-process NRMSE（missing=1） | complete-contribution NRMSE（missing=1） |
|---|---|---:|---:|---:|
| base_01 seed 6 | exact | 4/4 | 0.255 | 0.237 |
| sparse_01 seed 6 | 缺 B1 | 3/4 | 0.540 | 0.505 |
| sparse_02 seed 7 | exact | 4/4 | 0.688 | 0.311 |

其中 sparse_02 的 A/B1/B2 final loading norm ratio 分别为 7.009、3.669、0.932；对应四个 study-role 的 factor-process NRMSE 为 0.867、0.756、0.876、0.253，而完整贡献 NRMSE 为 0.320、0.287、0.397、0.238。这再次说明“载荷方向正确”不等于“过程尺度正确”。

全部 24 个端点的 study-role 覆盖为 67/96（0.698）；3 个 winner 为 11/12（0.917）；5 个 exact 端点为 20/20。Exact 端点的平均 factor-process NRMSE 为 0.340，平均完整贡献 NRMSE 为 0.252，明显优于全部端点的 0.594 和 0.502，但该结论仅来自 3 个自适应选择的数据集。

## 5. fixed-400 严格 practical convergence 的真正阻断项

24 个端点中，14 个 practical-converged，10 个为 slow case。第 400 步各门失败数如下：

| 门 | 全部失败 | 10 个 slow 中失败 |
|---|---:|---:|
| objective rate | 0 | 0 |
| objective monotonicity | 0 | 0 |
| short fitted | 7 | 7 |
| short RSS | 9 | 9 |
| short PPI quantile | 0 | 0 |
| short factor-PPI | 0 | 0 |
| long fitted | 9 | 9 |
| long RSS | 7 | 7 |
| long PPI quantile | 0 | 0 |
| long factor-PPI | 0 | 0 |

所以当前 fixed-400 批次的未收敛不能再归因于 PPI 门：PPI 门在 24/24 终点都通过。它也不是纯“语义误报”，因为多数 slow 端点存在可测的 fitted/RSS 移动；但移动幅度差异很大。例如 exact baseline seed 4 只以 short RSS `0.001047` 略高于 `0.001` 而失败，属于边界型；sparse_02 winner 的 long RSS 为 `0.03587`，属于真实慢移动。

5 个 exact 端点中只有 1 个 practical-converged，4 个为 slow；19 个非 exact 端点中 13 个 converged、6 个 slow。因此本批证据再次证明：convergence 状态不能替代结构正确性或科学质量判断。

## 6. 决策

Stage 4F 支持启动且只启动两条独立、版本化的 400→800 audit continuation：

1. sparse_02 seed 7：检验真实晚期尺度与 fitted/RSS 移动是否稳定或继续恶化；
2. base_01 seed 4：作为只略过 short-RSS 门的 exact baseline 对照。

不延长 sparse_01 winner；不把累计 800 结果用于替换 Stage 4E winner；不据此把全局主预算改为 800。后续方法优化应回到原有 6–9 时点场景，优先研究退火/早期普通 sweep 内的块竞争保护，以及 sparse dense-function 的尺度与报告映射，而不是先扩大大规模实验。
