# G12 靶向延长阶段 4G 完整报告

报告日期：2026-08-28
性质：开发性、真值后置解封、固定 horizon 的 400→800 continuation 机制审计
正式 v0L-V 结果：否
是否替换 Stage 4E 的 sweep-400 winner：否

## 1. 研究问题与边界

Stage 4G 只延长 Stage 4F 预先指定的两个 sweep-400 状态：

1. `g12ss4c_base_01__G12__04`：结构和身份完全正确的基线对照；在 sweep 400 仅 short-window RSS 略高于门槛。
2. `g12ss4c_sparse_02__G12__07`：Stage 4E 在真值隔离下按最大有效普通 T=1 ELBO 选出的 sparse_02 winner；结构和载荷方向正确，但 dense-grid、尺度和 fitted/RSS 仍存在明显问题。

每条轨迹从完整 sweep-400 变分状态继续 400 个普通 T=1 sweep，到累计 sweep 800。没有重新初始化、预得分、退火、winner 重选或真值参与。为观测完整固定 horizon，运行时把 `consecutive` 设为 401；因此运行器必然不会提前标记 practical convergence，最终是否满足公共 G12 的连续 5 步规则由离线审计判断。

本审计只有两个自适应选择的端点，回答的是机制问题，不估计推广成功率，也不能据此修改正式 v0L-V 或论文级成功率。

## 2. 执行、隔离与完整性

- continuation：2/2 COMPLETE；两条均 objective eligible、ELBO bridge 单调、源状态逐值一致性误差为 0；
- truth-free continuation：2/2 的 `truth_used=FALSE`；
- 两条 continuation 全部完成后才运行真值评价，wrapper 记录 `truth_read_after_continuations_complete=TRUE`；
- 评价：2 条 continuation、4 个 source/800 端点、16 个 study-role 真值行；QC 9/9 PASS；
- 没有重开 Stage 4E winner 选择，没有自动进一步延长，结果不属于正式 v0L-V；
- 两条各产生 21 个方向/尺度 checkpoint、完整 fit RDS、ELBO、practical diagnostics、warning、runtime 和 COMPLETE 标记；
- continuation 并发墙钟时间 1:42:08。基线耗时 6120.18 秒、峰值约 681.7 MiB；sparse_02 耗时 5326.94 秒、峰值约 401.7 MiB；总驱动器最大 RSS 约 790.7 MiB；
- 每条出现 1 个预期 warning：由于固定 horizon 使用 `consecutive=401`，400 个 continuation sweep 内不允许提前终止。这是协议性提示，不是拟合错误；
- 真值评价耗时 11.83 秒，退出码均为 0。

完整本地归档为 `v0lv-stage4g-complete-20260828-v1.tar.gz`，大小约 179 MiB；下载后 SHA256 与服务器一致：

`9fe85a72eff2c71eef99af3683e0af784cee9e6e9ea8aa230b98f64b962f1f45`

## 3. 数值收敛：两条轨迹到 sweep 800 均满足公共 G12 门槛

| 数据/端点 | sweep | final gate | last 5 全通过 | 公共 G12 反事实判定 | short fitted | short RSS | long fitted | long RSS |
|---|---:|---|---|---|---:|---:|---:|---:|
| base_01 seed 4 | 400 | FAIL | FAIL | FAIL | 0.000203 | 0.001047 | 0.001006 | 0.005498 |
| base_01 seed 4 | 800 | PASS | PASS | PASS | 0.000022 | 0.000122 | 0.000070 | 0.000382 |
| sparse_02 seed 7 | 400 | FAIL | FAIL | FAIL | 0.000502 | 0.003619 | 0.005111 | 0.035870 |
| sparse_02 seed 7 | 800 | PASS | PASS | PASS | 0.000036 | 0.000207 | 0.000115 | 0.000655 |

两条在 sweep 796–800 连续 5 步满足公共规则。objective、单调性、PPI quantile 和 factor-PPI 门槛均通过；400 时的阻断来自 fitted/RSS，800 时这些增量已经足够小。

这说明两个具体状态的“慢移动”最终会衰减，而不是无限振荡或 objective 失效。但它只证明数值增量变小，不能证明端点更接近真值。

## 4. 结构与身份：400→800 完全不变

两个端点在 sweep 400 和 800 均为：

- 选择因子数 `1,1,1`；
- 全局方向 3/3 correct；study-role 4/4 correct；
- missing、misplaced、extra、duplicate 均为 0；
- 8/8 feature components、4/4 covariance kernels 完成 loading-first matching；
- PPI 阈值 0.6、0.7、0.8、0.9 下结论完全一致；
- A、B1、B2 的 factor PPI 始终约为 0.998–1.000；载荷方向余弦始终高于 0.9991。

因此，增加 400 sweep 没有修复或破坏共享—特异身份。sparse_02 的问题不是末期发生了因子置换、符号错误或载荷方向错配，而是正确身份盆地内部的函数尺度、过程尺度和 dense-grid 外推质量问题。

## 5. 科学终点：ELBO 上升，但总体没有精度收益

### 5.1 基线 exact control

| 指标 | sweep 400 | sweep 800 | 800−400 / 相对变化 |
|---|---:|---:|---:|
| 普通 T=1 ELBO | -123663.752 | -123603.606 | +60.146 |
| 观测点 signal NRMSE | 0.186251 | 0.186187 | -0.000065（-0.035%） |
| dense-grid signal NRMSE | 0.321415 | 0.329533 | +0.008117（+2.53%） |
| loading relative L2 error | 0.146972 | 0.156693 | +0.009721（+6.61%） |
| matched function recall | 0.869604 | 0.867208 | -0.002397 |
| trajectory absolute correlation | 0.971767 | 0.969834 | -0.001933 |
| FPCA score absolute correlation | 0.955732 | 0.956033 | +0.000301 |
| feature total ISE | 0.171159 | 0.175043 | +2.27% |
| estimate-to-projection ISE | 0.102313 | 0.107121 | +4.70% |
| kernel relative ISE | 0.163608 | 0.165871 | +1.38% |
| covariance-operator relative error | 0.519228 | 0.537647 | +3.55% |

基线的观测点重构只改善约 0.035%，其余多数科学指标轻微变差。800 端点的 dense-grid NRMSE 仍明显优于 mean-only 0.674，但没有优于 400 端点。

### 5.2 sparse_02 exact winner

| 指标 | sweep 400 | sweep 800 | 800−400 / 相对变化 |
|---|---:|---:|---:|
| 普通 T=1 ELBO | -93033.062 | -93004.139 | +28.923 |
| 观测点 signal NRMSE | 0.211163 | 0.211231 | +0.000069（+0.032%） |
| dense-grid signal NRMSE | 3.579743 | 3.661408 | +0.081664（+2.28%） |
| loading relative L2 error | 2.919554 | 3.006811 | +0.087257（+2.99%） |
| matched function recall | 0.503770 | 0.503170 | -0.000600 |
| trajectory absolute correlation | 0.954336 | 0.953063 | -0.001273 |
| FPCA score absolute correlation | 0.506375 | 0.503176 | -0.003198 |
| feature total ISE | 0.896051 | 0.904291 | +0.92% |
| estimate-to-projection ISE | 0.772004 | 0.781338 | +1.21% |
| kernel relative ISE | 1.478659 | 1.502866 | +1.64% |
| covariance-operator relative error | 40.134599 | 42.636373 | +6.23% |

sparse_02 的观测点误差几乎不变，但 dense-grid、载荷、feature、score、kernel 和协方差指标全部没有改善。800 端点 dense-grid NRMSE 3.661，远差于同端点 mean-only 0.706。普通 ELBO 的继续增加与未观测时间区域的科学恢复方向相反。

## 6. 过程、完整贡献与尺度轨迹

### 6.1 study-role 汇总

| 数据 | sweep | matched/total | factor-process NRMSE | complete-contribution NRMSE | mean canonical loading-norm ratio |
|---|---:|---:|---:|---:|---:|
| base_01 | 400 | 4/4 | 0.254411 | 0.234773 | 1.013540 |
| base_01 | 800 | 4/4 | 0.261048 | 0.243162 | 1.023357 |
| sparse_02 | 400 | 4/4 | 0.687910 | 0.310523 | 4.654785 |
| sparse_02 | 800 | 4/4 | 0.689597 | 0.331140 | 4.766296 |

基线 process NRMSE 增加 2.61%，完整贡献 NRMSE 增加 3.57%；主要由研究 2 特异因子的 contribution NRMSE 从 0.2970 增至 0.3261 驱动。sparse_02 process NRMSE 只增加 0.25%，但完整贡献 NRMSE 增加 6.64%；其中研究 1 的共享贡献误差从 0.3205 增至 0.4034，是主要恶化来源。

### 6.2 底层尺度继续移动，但方向保持不变

基线原始 CAVI loading-norm ratio 从 A/B1/B2=`1.584/1.202/1.102` 移到 `1.143/0.974/0.923`，更接近 1；同时 factor-process grid RMS 从 `0.602/0.748/1.188` 增到 `0.835/0.925/1.460`。complete-contribution RMS proxy 只从 `0.286/0.270/0.393` 变为 `0.287/0.270/0.404`。这主要表现为载荷—过程之间持续但逐渐减慢的尺度重分配，而非结构改变。

sparse_02 的原始 loading-norm ratio 也从 `1.457/1.100/1.217` 移到 `1.200/0.950/1.042`，方向余弦仍约为 0.999；但 A/B1/B2 的 feature grid RMS 从 `3.162/2.654/0.651` 增至 `4.596/2.867/0.720`，factor-process grid RMS 从 `5.465/3.644/0.780` 增至 `6.817/4.313/0.910`。canonical factor scale 进一步增加：

- A：4.811→5.996（+24.6%）；
- B1：3.335→3.942（+18.2%）；
- B2：0.766→0.893（+16.6%）。

因此 sparse_02 的晚期移动不是寻找另一结构盆地，而是当前正确身份盆地内部继续扩大时间函数/因子过程尺度。载荷与得分的逆向尺度补偿使观测点拟合和部分贡献 proxy 看似稳定，却不能约束未观测 dense-grid 上的函数幅度与协方差。

## 7. 结论

1. **800 sweep 可以让这两条特定轨迹通过数值 practical gate，但不能提升科学恢复。** 两条在 sweep 796–800 连续 5 步通过公共 G12 门槛，说明 fitted/RSS 慢移动最终衰减。
2. **不应把全局最大预算改成 800。** 额外 400 sweep 约增加一倍计算成本；结构完全不变，观测重构几乎不变，dense-grid、载荷、feature、covariance 和贡献总体略有或明显恶化。
3. **严格 practical convergence 不是正确性指标。** sparse_02 最终在错误很大的 dense-grid/协方差尺度上数值稳定；“收敛”只表示相邻迭代变化小，不表示端点接近真值。
4. **sparse_02 的核心是弱时间信息下的函数尺度/外推失控，不是因子身份或 PPI。** PPI、因子数、载荷方向和 matching 全部正确且稳定，继续迭代没有针对性。
5. **Stage 4E 的 sweep-400 winner 和当前 G12 预算不变。** Stage 4G 是后验机制证据，不能替换原 winner、重估成功率或进入正式 v0L-V。

## 8. 下一步建议

### 8.1 立即固化，不再延长

- 把本报告、两条完整 continuation、评价表、QC、日志和归档 SHA 作为 Stage 4G v1 证据保留；
- 不运行 800→更长，不把 800 设为默认预算；
- G12 仍作为下一版开发协议的首选初始化；历史 `current_1_over_m` 接口与既有冻结实验保持不变。

### 8.2 阶段 5A：回到原始 6–9 时点场景，定位早期盆地分流

优先回答 Stage 4F 尚未定位的问题：失败块究竟在 99 个退火 sweep 的哪一步、哪一参数更新后被关闭。只增加轻量 checkpoint，不改变模型或 ELBO；记录每个退火 sweep及最初 20 个普通 sweep 的 block PPI、loading energy、score/process energy、RSS/ELBO 和共享—特异重叠。使用少量已有原始密度数据和预注册拟合种子，不先扩大 Monte Carlo。

### 8.3 阶段 5B：只比较一个可解释、低成本的优化日程候选

若 5A 证明确为“连续参数尚未形成就被 spike-slab 关闭”，再比较当前 G12 与单一候选，例如短暂的 selection burn-in：先允许载荷、feature 和 score 建立表示，再开启 PPI 竞争。它只改变优化日程，不改变模型、先验、普通 ELBO 或真值盲 winner 选择；必须在相同起点数和相同总 sweep 预算下比较，避免混合多种技巧。

判定门槛应预先限定为：真值盲 winner 的结构命中率、独立起点支持、有效普通 ELBO、观测/dense 重构、载荷—feature—process—贡献链和计算成本。若没有跨数据集增益，保留现有 G12。

### 8.4 中期问题

- 在原始 6–9 时点上审计 dense-grid 与观测点误差是否仍分离；若只在 3–5 时点极端稀疏场景出现，则把它界定为适用边界，而不是立刻改模型；
- 若原始密度下也出现尺度失控，再研究时间平滑/尺度识别约束与后处理映射，但这属于后续模型版本，不能从本两端点直接推出；
- 单独设计跨 `p` 和因子强度的维数选择实验，避免仅用高维 union factor-PPI 判断因子存在；
- 核心无协变量路径稳定后，再验证 `d>0` 的协变量效应；当前证据不能外推到协变量版本。

## 9. 关键机器结果

- `ALL_2_TERMINALS.csv`
- `evaluation_20260828_v1/EVALUATION_QC.csv`
- `evaluation_20260828_v1/SOURCE_AND_800_GATE_AUDIT.csv`
- `evaluation_20260828_v1/SOURCE_AND_800_ENDPOINT_METRICS.csv`
- `evaluation_20260828_v1/SOURCE_TO_800_METRIC_CHANGES.csv`
- `evaluation_20260828_v1/CONTINUATION_MECHANISM_TIMELINE.csv`
- `evaluation_20260828_v1/SOURCE_AND_800_PROCESS_CONTRIBUTION.csv`
- `evaluation_20260828_v1/SOURCE_AND_800_PROCESS_CONTRIBUTION_SUMMARY.csv`
