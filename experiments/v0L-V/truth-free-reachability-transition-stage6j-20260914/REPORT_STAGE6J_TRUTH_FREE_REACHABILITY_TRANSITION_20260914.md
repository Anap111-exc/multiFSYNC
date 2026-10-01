# Stage 6J：真值盲可达性转折点与 residual-FPCA 配对审计报告

报告日期：2026-09-14  
协议：`STAGE6J_TRUTH_FREE_REACHABILITY_TRANSITION_V1_20260914`  
性质：开发期、单个已开发困难数据种子的机制诊断；不是论文级 Monte Carlo，也不是跨维度选择验证。

## 1. 结论摘要

Stage 6J 在同一固定真实维度、同一数据和 4 个配对拟合种子下比较了：

1. 当前首选的随机 Gram 单位函数能量初始化（历史简称 G12）；
2. 真值盲、study-aware residual-FPCA 初始化。

两臂均保持一次预得分、Jaoua 退火和 400 个普通 `T=1` sweep 不变。8/8 个端点完成、有效且通过 practical convergence，0 error、0 warning。真值在端点及最大普通 ELBO winner 冻结后才解封。

主要结果为：

- G12 core 恢复 `0/4`，residual-FPCA core 恢复 `0/4`；
- residual-FPCA 的 4 个端点全部塌缩为 `(shared,B1,B2)=(0,0,0)`；
- G12 的 3 个端点为 `(1,0,0)`，另 1 个为身份错置的 `(1,0,1)`；
- 最大普通 ELBO 的真值盲 winner 是 `s6j__G12__02`，但它只恢复共享方向，两个研究特异方向均 missing；
- residual-FPCA 的数据驱动初始化确实执行成功、无 fallback，但三个载荷块在最初 2--5 个退火 sweep 内全部归零；失败发生在早期优化路径，而非终点选择；
- G12 在 ELBO、重构、协方差和科学结构覆盖上均明显优于 residual-FPCA，因此 residual-FPCA 不能替换 G12，也不支持再叠加为混合默认路径。

结合 Stage 6H 的同数据 12 个完整数据 G12 起点 `0/12` core，以及 Stage 6I 的 6 条 exact/near-truth 轨迹 `6/6` core，可以作出更准确的定位：

> 当前模型和普通 ELBO在该数据上包含稳定且目标值更高的正确盆地，但现有两种真值盲初始化均未生成正确端点；当前首要阻断仍是候选盆地可达性，不是新留出准则造成，也不是普通 ELBO在已经出现的正确与错误端点间选错。

以上结论只直接适用于数据种子 `84611001`。它是机制反例，不能估计 G12 在一般数据种子上的失败率。

## 2. 设计、数据与隔离

- 数据：只读复用 Stage 6G 的 `s6g_minimal_L111_M2` 观测包；数据种子 `84611001`；
- `S=2`、`n_s=(30,30)`、`p=500`、`d=0`、每人 6--9 个不规则观测时点；
- 真值与拟合维度：`L_f=1`、`L_s=(1,1)`，所有因子的 `M=2`，`K=5`；
- 拟合种子：`84611101:84611104`，两臂严格配对；
- 共同路径：一次预得分、Jaoua 退火 `c(1,1.9,100)`、99 个退火 sweep、400 个普通 `T=1` sweep；
- `lambda_orth=0`，每条 fit 内 `n_cpus=1`，8 条独立 fit 在任务层并发；
- 无 racing、无 continuation、无跨维度候选、无 one-SE 或预测准则选择；
- 唯一真值盲 endpoint rule：8 个有效端点中最大普通 `T=1` ELBO；
- truth 未用于拟合、停止、eligibility 或 winner 选择。

Stage 6J 不是重复证明“G12 在其他数据上能找到正确端点”，而是针对一个已由 Stage 6H 确认 `12/12` G12 均失败、又由 Stage 6I oracle 证明正确盆地存在的困难数据，检验 residual-FPCA 是否能改变候选生成机制。

## 3. 运行与完整性

| 项目 | 结果 |
|---|---:|
| 完成 fit | 8/8 |
| objective eligible | 8/8 |
| practical convergence | 8/8 |
| fit error | 0 |
| fit warning | 0 |
| 科学评价成功 | 8/8 |
| 科学评价 error | 0 |
| continuation | 0 |
| 新数据生成 | 0 |
| 跨维度选择 | 未运行 |

8 条 fit 的累计单核时间约 4.32 小时；8 路并发墙钟约 34.1 分钟。G12 单条中位时间约 33.2 分钟，residual-FPCA 约 31.5 分钟。residual-FPCA 约快 5.6%，但这一小幅速度差不能补偿其结构完全塌缩。

## 4. 真值盲终点与 ELBO

| 初始化 | seed | 普通 ELBO | 保留数 `(A,B1,B2)` | core | 观测 NRMSE | dense NRMSE | covariance error |
|---|---:|---:|---|---:|---:|---:|---:|
| G12 | 1 | -114970.1060 | `(1,0,0)` | 否 | 0.16051 | 0.18487 | 0.36481 |
| G12 | 2 | **-114962.7461** | `(1,0,0)` | 否 | 0.16051 | 0.18493 | 0.36636 |
| G12 | 3 | -114982.3680 | `(1,0,0)` | 否 | 0.16051 | 0.18369 | 0.40550 |
| G12 | 4 | -116425.0038 | `(1,0,1)` | 否 | 0.18002 | 0.20061 | 0.73779 |
| residual-FPCA | 1 | -118140.8141 | `(0,0,0)` | 否 | 0.21491 | 0.23287 | 1.00000 |
| residual-FPCA | 2 | -118140.8141 | `(0,0,0)` | 否 | 0.21491 | 0.23287 | 1.00000 |
| residual-FPCA | 3 | -118140.8141 | `(0,0,0)` | 否 | 0.21491 | 0.23287 | 1.00000 |
| residual-FPCA | 4 | -118140.8141 | `(0,0,0)` | 否 | 0.21491 | 0.23287 | 1.00000 |

4 个 residual-FPCA 端点具有不同 fit 种子、不同 independence ID 和不同 fit SHA256；初始化诊断也有轻微差异。因此相同的终点数值不是同一文件被重复计数，而是 4 条轨迹进入了同一个强吸引的空因子盆地。

G12 对 residual-FPCA 的配对优势在 4/4 个种子上方向一致：ELBO 高 `1715.8--3178.1`，观测和 dense NRMSE 均更低。按臂均值，G12 相比 residual-FPCA：

- 观测点 signal NRMSE 低约 23.0%；
- dense-grid signal NRMSE 低约 19.0%；
- covariance operator 相对误差低约 53.1%；
- matched feature recall 为 0.941，而 residual-FPCA 因没有任何匹配成分为 0。

residual-FPCA 的 loading error 是 `NA`，含义是没有可匹配的载荷方向，不能解释为零误差。

## 5. 真值盲 winner 的科学质量

winner `s6j__G12__02` 的结构为 `(1,0,0)`：

- shared factor PPI `R_A=0.99229`；
- 两个 specific factor PPI 均为 0；
- 真实方向：1 个 correct、2 个 missing；
- study-role：2/4 匹配；
- FPCA 特征成分：4/8 匹配；
- exact structure、core 和完整特征覆盖均失败。

已恢复共享方向的 conditional 科学质量较好：

| 指标 | winner |
|---|---:|
| 观测点 signal NRMSE | 0.16051 |
| dense-grid signal NRMSE | 0.18493 |
| dense relative MISE | 0.03420 |
| covariance operator 相对误差 | 0.36636 |
| conditional loading relative L2 | 0.10301 |
| loading-PPI AUC | 0.95716 |
| matched function recall | 0.95565 |
| matched trajectory abs(cor) | 0.92741 |
| matched score abs(cor) | 0.91225 |
| matched kernel relative ISE | 0.04748 |

但把 missing 方向按零估计计入后：

- factor-process NRMSE 为 0.69403；
- complete-contribution NRMSE 为 0.70268。

因此，该 winner 是“共享方向质量很好但共享—特异结构不完整”的解。较好重构、较高 ELBO、practical convergence 和低 conditional error 均不能替代完整身份与覆盖评价。

G12 seed 4 的 `(1,0,1)` 也不是部分正确的 `shared + B2`：两个保留方向发生共享/特异身份错置，结果为 0 correct、2 misplaced、1 missing。

## 6. residual-FPCA 为什么失败

初始化本身在技术上完成：

- 4/4 条均报告 `used_method=residual_fpca`；
- 3 个候选块均成功，无数值秩 fallback；
- 总解释比例约 0.22182--0.22191；
- 初始化耗时约 2.1 秒。

然而，dense driver trace 显示它没有在 CAVI 中存活：

| 块 | sweep 1 posterior norm | sweep 2 posterior norm | sweep 5 posterior norm |
|---|---:|---:|---:|
| shared loading | 0.1661 | 0.0000064 | 0 |
| study 1 specific loading | 0.8277 | 0.2357 | 0 |
| study 2 specific loading | 0.3701 | 0.00334 | 0 |

第 5 个 sweep 时，相关载荷、得分和时间函数的数据 driver 已全部为零。之后 PPI 在退火先验动态中短暂保持中间值，但到 sweep 99 已降至约 `10^-10`，三个块此后不可逆地停留在空盆地。换言之：

1. 这不是 400 sweep 不够；空盆地在前 5 个退火 sweep 已形成；
2. 这不是最终 ELBO selector 删除了有用候选；有用候选在形成终点前已经消失；
3. 单独构造更“数据驱动”的均值初始化不够，初始化尺度、得分/时间函数更新与早期稀疏反馈必须相容；
4. residual-FPCA 不应成为下一版默认初始化，也没有证据支持把它与 G12 叠加为昂贵且难解释的混合方案。

现有 trace 能定位到早期退火/CAVI 闭环，但还不能仅凭终点断言是预得分、载荷更新尺度、`Beta(1,p)` 反馈或退火温度中的哪一个单独造成。任何进一步试验都应是一次小型、配对、带停止规则的时序/稀疏反馈因果审计，而不是继续堆初始化器。

## 7. G12 的失败路径

G12 没有整体塌为空盆地，而是稳定进入 shared-only 盆地：

- shared loading 在 sweep 1 的平均 posterior norm 为 19.62，factor PPI 很快达到 1 并持续存活；
- study 1 specific 从 sweep 1 的 6.17 降至 sweep 5 的 0.277；sweep 75 后快速关闭，sweep 99 归零；
- study 2 specific 的证据维持更久，sweep 50--75 仍有较强局部 PPI，但随后逐渐衰减，到总 sweep 499 平均 PPI 约 0.0199，未达到结构保留；
- 这与 Stage 6H 中物理 `B10_1` 方向 36/36 missing 的方向性失败一致。

因此 G12 对该数据的故障不是初始化瞬间把所有块设为零，而是共享块在顺序更新和稀疏反馈中占优，两个特异块在退火期间逐步失去证据。

## 8. 与 Stage 6H/6I 的联合解释

三阶段共同形成以下证据链：

| 阶段 | 干预 | 结果 | 能回答的问题 |
|---|---|---|---|
| 6H | G12 起点增至 12；恢复完整 6--9 时点 | 完整数据仍 0/12 core | 失败不能仅归因于 4 起点或留出数据损失 |
| 6I | exact/near-truth oracle；研究顺序交换 | oracle 6/6 core，正确盆地 ELBO更高 | 模型/普通目标中存在稳定正确盆地；非确定性 study-1 位置 bug |
| 6J | G12 对 residual-FPCA，4 个配对种子，密集早期 trace | 两臂均 0/4；residual 在 sweep 5 前全塌缩 | residual-FPCA 不能解决可达性；阻断发生在早期路径 |

Stage 6I 正确盆地最大普通 ELBO约为 `-114702.3826`：

- 比 Stage 6J 四起点 G12 winner 高约 260.36；
- 比 Stage 6H 十二起点 G12 winner（`-114892.9955`）高约 190.61；
- 比 residual-FPCA 空盆地高约 3438.43。

所以，在同一真实固定维度内，若真值盲路径生成了 Stage 6I 所证明的正确端点，普通 ELBO应能优先于现有错误端点。当前没有证据要求替换固定维度内的最大普通 ELBO规则。

这不代表过指定 `L/M` 时的跨维度问题已经解决：不同维度间的普通 ELBO不可直接承担同一含义，留出风险/joint log-score/one-SE 的候选准则仍需独立验证。但 Stage 6G 在这个数据上的判别失败不能单独归因于准则，因为每个候选配置首先缺少正确盆地。

## 9. 本阶段可以与不能得出的结论

### 可以得出

1. residual-FPCA 在该困难数据上明显差于 G12，不能晋升为默认初始化。
2. residual-FPCA 的失败是前 2--5 个退火 sweep 内的动态塌缩，不是最终 selector 或迭代不足。
3. G12 对该数据的失败主要是 shared-only 竞争盆地；增加到 12 个起点在 Stage 6H 已经验证仍不能恢复 core。
4. 当前困难数据的首要问题是正确候选生成/可达性；新维度准则没有制造这一失败。
5. practical convergence 只表示数值停止，不表示科学结构正确。

### 不能得出

1. 不能由一个开发数据种子宣称 G12 在一般场景失效。
2. 不能把 residual-FPCA 的失败推广为所有 PCA/谱初始化必然无效。
3. 不能仅凭本轮断言 `Beta(1,p)`、退火、预得分、更新顺序或某个先验是唯一根因。
4. 不能用本轮评价跨维度预测准则，因为本轮没有跨维度候选。
5. 不能把 oracle 初始化用于正式拟合或方法选择；它只证明正确盆地存在。

## 10. 建议的后续方向与停止规则

优先级应保持克制：

1. 固化 Stage 6J 证据，继续保留 G12 + 一次预得分 + Jaoua 退火 + 多起点；不采用 residual-FPCA 或混合初始化。
2. 若继续追查可达性，只允许一次预注册、很小规模的早期稀疏反馈/退火时序配对实验。它应直接检验“初始数据证据是否因过早 gate/PPI 反馈而被不可逆关闭”，而不是增加更多通用初始化方案。
3. 若该最小干预不能在同一困难数据上产生 core，则把这一数据种子记录为现有算法局限，停止继续在它上面追加补丁。
4. 论文主线转回导师要求的两项工作：明确的 `L/M` 外部选择准则模拟，以及固定正确维度下的协变量效应模拟。
5. 维度准则实验必须先报告候选池可达性；没有正确端点时应记为“candidate-generation failure”，不能把责任归给跨维度 selector。

## 11. 结果落盘与服务器状态

服务器完整结果保留在：

`/root/multifsync-results/stage6j_20260914_v1`

本地版本化证据位于本报告同目录：

- `stage6j_server_compact_20260914_v1/`：CSV/文本紧凑证据；
- `stage6j_server_compact_20260914_v1.tar.gz`；
- 归档 SHA256：`42F2D6E450B16E8CB76D6804603713BD889842262620BB1311565F8898BEB7F7`；
- `STAGE6J_TRUTH_UNSEAL_AUTHORIZATION_20260914.txt`；
- `evaluate_stage6j_20260914_v1.R`。

评价只读取已冻结的 8 个端点，没有启动新 fit、continuation、数据生成或跨维度选择。结果下载和完整性核验后，ECS 的 `StopCharging` 请求已被 API 接受。
