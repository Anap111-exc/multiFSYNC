# G12 正确盆地科学质量只读审计 Stage 5A 完整报告

报告日期：2026-08-29
性质：开发性、真值后置、只读审计
正式 v0L-V 结果：否
是否运行新数据、拟合或 continuation：否

## 1. 结论摘要

1. 本次审计完整继承两个历史实验中原始 6--9 个不规则时点的 104 条 G12 端点：10 个数据集、45 条正确共享--特异身份端点、10 个冻结的真值盲最大普通 ELBO winner。输入和 QC 全部通过。
2. 单起点正确身份命中仍有明显数据间和起点间波动：45/104（43.27%）；各数据为 25%--75%。因此 G12 改善了可达性，但没有消除非凸多盆地，多起点仍不可删除。
3. 10/10 个冻结 winner 均为正确身份；每个数据的最佳正确端点 ELBO 均高于最佳错误端点，差值为 4561.9--9014.4（中位数 6807.6）。这为“有效普通 ELBO + 多起点”提供了稳定的真值盲选择支持。
4. 正确身份并不自动保证相同的科学质量。45 条正确端点中存在少量尺度异常子盆地，最严重端点的 dense-grid signal NRMSE、loading relative L2 error 和 covariance-operator relative error 分别达到 3.221、2.078 和 30.780；但对应的冻结 winner 最大值降为 0.531、0.353 和 1.234，灾难性尾部未进入 winner。
5. 异常端点并非载荷方向、因子身份或 component coverage 错误：45/45 条正确端点均完成 8/8 feature component 和 4/4 kernel 匹配；主队列 27 条正确端点的 108/108 个 study-role 过程也均按载荷身份匹配。主要异常是方向近乎正确时的尺度膨胀及其 dense-grid、协方差、过程和贡献后果。
6. 最大 ELBO 对观测点重构具有较强的一致性，但不是所有真值指标的 oracle。正确盆地内，ELBO 与观测点 NRMSE 的数据内 Spearman 中位数为 -0.90，9 个可估数据中 8 个方向一致；对 dense、loading、feature、kernel 和 covariance 的一致性较弱。保留 ELBO 作为真值盲选择准则是合理的，但不能宣称它逐项优化科学恢复。
7. 104/104 条端点均运行到 200 sweep，0/104 达到 strict practical convergence；因此严格收敛在本批结果中没有可比较变异，不能充当结构或科学成功指标。辅助 endpoint pass 在正确身份端点内对应更好的科学指标中位数，但它只覆盖 8/45 条正确端点，同时 29 条 pass 中有 21 条身份不正确，不能替代 ELBO、身份评价或多起点。
8. Stage 5A 不给出修改模型、先验、ELBO、G12、退火或预算的证据。当前最合适的解释仍是：G12 + 预得分 + Jaoua 退火 + 12 起点 + 最大有效普通 ELBO，是已经得到支持的优化缓解方案；优化路径稳定性得到管理，但尚未被理论或算法上彻底解决。

## 2. 问题、输入与证据边界

Stage 5A 不重复 G12 拟合，而回答四个只读问题：

1. 已进入正确共享--特异身份盆地的端点，其科学恢复是否同质；
2. 冻结的真值盲最大 ELBO winner 是否能够避开科学质量差的正确身份子盆地；
3. 原始 6--9 时点场景是否仍出现尺度、dense-grid、过程或完整贡献异常；
4. strict practical convergence 是否具有可用的结构或科学质量含义。

输入分为两个保持独立身份的历史队列：

- 主队列：`g12ic_01`--`g12ic_06`，数据种子 `82622001`--`82622006`，强可辨识，12 个 G12 起点，共 72 条；本地具有完整 fit/evaluation/truth RDS，可重算过程、完整贡献和尺度。
- 次队列：`g12cs_base_01/02` 与 `g12cs_weak_01/02`，数据种子 `82623001`--`82623004`，每数据 8 个 G12 起点，共 32 条；本地只有冻结汇总结果，不能重算过程和完整贡献。

104 条起点不是 104 个独立数据重复。主要分析单位是数据集；起点级统计用于描述同一数据内的盆地异质性。两个历史队列也没有被重新包装成一个预注册推断实验，合并数字只用于开发性描述。

## 3. 执行完整性

- 输入数据集：10；G12 端点：104；fit_id：104 个唯一值。
- 主队列完整 RDS：72/72；次队列冻结汇总：32/32。
- 正确身份端点：45；冻结 winner：10；正确 winner：10。
- objective eligible：104/104；`truth_used_for_fit_or_selection=FALSE`：104/104。
- 主队列重算：288 个 study-role 过程/贡献行、216 个 loading-scale block 行。
- QC：15/15 PASS；完成标记为 `status=PASS`。
- 新数据生成、新拟合、continuation、winner 重选：全部为 FALSE。

## 4. 正确盆地可达性与 ELBO 支持

### 4.1 单起点结构结果

104 条端点可分为：

- 正确因子数且正确身份：45（43.27%）；
- 因子数为 1/1/1、但共享--特异身份错误：30（28.85%）；
- 缺失因子或其他不完整结构：29（27.88%），其中 1/0/0 为 2 条、1/0/1 为 12 条、1/1/0 为 15 条。

因此，75/104 条的因子数表面上是 1/1/1，但其中 30 条身份错误。因子数或 factor-PPI 不能代替载荷优先的一对一身份匹配。

主队列正确率为 27/72（37.50%）；次队列为 18/32（56.25%）。若仅作描述，8 个 baseline/strong 数据合计 34/88（38.64%），两个 weak-separation 数据合计 11/16（68.75%）。后者只有两个数据集且来自既有跨场景实验，不能据此宣称弱可辨识场景更容易。

### 4.2 数据级结果

| 数据 | 场景 | 起点 | 正确身份 | 数量对但身份错 | 不完整 | 最佳正确减最佳错误 ELBO |
|---|---|---:|---:|---:|---:|---:|
| g12cs_base_01 | baseline | 8 | 2（25.0%） | 1 | 5 | 7065.0 |
| g12cs_base_02 | baseline | 8 | 5（62.5%） | 2 | 1 | 9014.4 |
| g12cs_weak_01 | weak separation | 8 | 6（75.0%） | 1 | 1 | 6185.8 |
| g12cs_weak_02 | weak separation | 8 | 5（62.5%） | 1 | 2 | 4561.9 |
| g12ic_01 | baseline | 12 | 7（58.3%） | 1 | 4 | 6702.3 |
| g12ic_02 | baseline | 12 | 4（33.3%） | 5 | 3 | 7570.8 |
| g12ic_03 | baseline | 12 | 3（25.0%） | 4 | 5 | 8159.1 |
| g12ic_04 | baseline | 12 | 3（25.0%） | 7 | 2 | 6871.4 |
| g12ic_05 | baseline | 12 | 5（41.7%） | 3 | 4 | 5619.7 |
| g12ic_06 | baseline | 12 | 5（41.7%） | 5 | 2 | 6743.7 |

10 个数据中，最佳正确端点的 ELBO 均严格高于最佳错误端点；冻结 winner 也逐一等于各数据的最大有效普通 ELBO 端点。这个结果支持当前选择规则，但不能把开发数据上的 10/10 直接解释为未来论文级数据上的成功概率。

## 5. 正确身份盆地内部的科学质量

### 5.1 45 条正确端点与 10 个冻结 winner

下表为端点级中位数 `[最小值, 最大值]`。过程和完整贡献只在具有完整 RDS 的主队列可得，因此相应分母为 27 条正确端点和 6 个 winner；其余指标分母为 45 和 10。

| 指标 | 全部正确身份端点 | 冻结 winner |
|---|---:|---:|
| 观测点 signal NRMSE | 0.168 `[0.136, 0.226]` | 0.167 `[0.136, 0.195]` |
| dense-grid signal NRMSE | 0.290 `[0.168, 3.221]` | 0.280 `[0.169, 0.531]` |
| loading relative L2 error | 0.138 `[0.097, 2.078]` | 0.138 `[0.106, 0.353]` |
| feature total ISE | 0.097 `[0.027, 0.568]` | 0.090 `[0.028, 0.332]` |
| feature estimate-to-projection ISE | 0.041 `[0.007, 0.517]` | 0.038 `[0.008, 0.272]` |
| kernel relative ISE | 0.191 `[0.038, 0.723]` | 0.187 `[0.039, 0.374]` |
| covariance-operator relative error | 0.613 `[0.342, 30.780]` | 0.551 `[0.346, 1.234]` |
| matched function recall mean | 0.922 `[0.625, 0.976]` | 0.923 `[0.790, 0.975]` |
| trajectory absolute correlation | 0.973 `[0.923, 0.989]` | 0.975 `[0.960, 0.989]` |
| FPCA score absolute correlation | 0.949 `[0.699, 0.990]` | 0.953 `[0.832, 0.989]` |
| factor-process NRMSE | 0.247 `[0.196, 0.605]` | 0.245 `[0.196, 0.293]` |
| complete-contribution NRMSE | 0.214 `[0.155, 0.344]` | 0.210 `[0.155, 0.269]` |

中位数层面，正确身份端点总体科学质量良好；真正需要注意的是 dense、loading、covariance 和 process 的右尾。冻结 winner 的中位数通常相近或略好，但更重要的是明显截短了这些灾难性尾部。

### 5.2 每个数据的冻结 winner

| 数据 | obs NRMSE | dense NRMSE | loading error | feature ISE | kernel ISE | covariance error | trajectory cor | score cor |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| g12cs_base_01 | 0.186 | 0.431 | 0.192 | 0.237 | 0.286 | 0.843 | 0.969 | 0.871 |
| g12cs_base_02 | 0.137 | 0.169 | 0.112 | 0.028 | 0.039 | 0.346 | 0.987 | 0.989 |
| g12cs_weak_01 | 0.163 | 0.400 | 0.251 | 0.119 | 0.180 | 0.940 | 0.973 | 0.949 |
| g12cs_weak_02 | 0.190 | 0.531 | 0.353 | 0.332 | 0.374 | 1.234 | 0.960 | 0.832 |
| g12ic_01 | 0.166 | 0.224 | 0.107 | 0.072 | 0.121 | 0.424 | 0.980 | 0.956 |
| g12ic_02 | 0.181 | 0.287 | 0.106 | 0.097 | 0.213 | 0.473 | 0.976 | 0.949 |
| g12ic_03 | 0.136 | 0.190 | 0.138 | 0.064 | 0.054 | 0.421 | 0.989 | 0.980 |
| g12ic_04 | 0.156 | 0.192 | 0.139 | 0.072 | 0.118 | 0.420 | 0.980 | 0.957 |
| g12ic_05 | 0.195 | 0.358 | 0.135 | 0.211 | 0.312 | 0.669 | 0.965 | 0.923 |
| g12ic_06 | 0.168 | 0.274 | 0.152 | 0.084 | 0.194 | 0.630 | 0.974 | 0.962 |

主队列 6 个 winner 的 factor-process NRMSE 为 0.196--0.293，complete-contribution NRMSE 为 0.155--0.269。次队列没有本地完整 RDS，因此没有补造这两项结果。

`g12cs_weak_02` 是 10 个 winner 中科学质量最弱的一条：dense NRMSE 0.531、loading error 0.353、covariance error 1.234。它仍远好于主队列中未被选择的灾难性正确端点，但表明“winner 避开极端失败”不等于“每个 winner 的所有指标都同样优秀”。

## 6. 异常机制：身份正确但尺度错误

45 条正确端点全部完成 8/8 feature components 和 4/4 covariance kernels 的一对一匹配；主队列 27 条正确端点的 108/108 个 study-role 真值项也全部按 loading identity 匹配。因此右尾不是 feature-missing、kernel-missing 或共享--特异错配造成的。

主队列载荷审计显示，正确端点的方向余弦始终非常高：A、B1、B2 的最小 canonical loading cosine 分别为 0.999717、0.999309 和 0.999484。与之相反，canonical loading norm ratio 存在明显右尾：

| block | 全部正确端点 norm ratio 中位数 / 最大值 | winner 中位数 / 最大值 | 全部正确端点 factor scale 中位数 / 最大值 |
|---|---:|---:|---:|
| A | 1.085 / 6.615 | 1.011 / 1.097 | 0.477 / 2.740 |
| B1 | 1.136 / 2.740 | 1.011 / 1.264 | 0.604 / 1.750 |
| B2 | 1.099 / 1.492 | 1.112 / 1.294 | 0.602 / 0.732 |

最典型端点为 `g12ic_06__G12IC__gram_unit_energy__06`：

- 共享 A 方向余弦为 0.999852，但 factor scale 为 2.740、canonical loading norm ratio 为 6.615；
- 观测点 NRMSE 仅升至 0.226，但 dense NRMSE 为 3.221；
- loading error 2.078、covariance error 30.780；
- factor-process NRMSE 0.605、complete-contribution NRMSE 0.344；
- 普通 ELBO 为 -127520.5，低于同数据冻结 winner 的 -123120.3，因此没有被选择。

同一数据的另外两个正确端点也出现较轻的 A 尺度膨胀：dense NRMSE 为 1.321 和 1.161，covariance error 为 7.233 和 5.660。冻结 winner 在多数科学指标上是该数据 5 个正确端点中的第 1 或第 2 名，避开了整组尺度右尾。

过程/贡献的 study-role 结果也支持同一解释。全部正确端点中，共享 A 在两个研究的 process NRMSE 最大值为 0.873 和 0.896，而冻结 winner 中相应最大值降为 0.345 和 0.329；B1/B2 winner 的均值约为 0.245/0.254。异常主要表现为正确方向上的尺度和时间函数幅度问题，而不是重新进入错误身份盆地。

## 7. 最大普通 ELBO 在正确盆地内选择了什么

在每个数据的正确端点内部，计算 ELBO 与科学指标的 Spearman 相关；至少需要 3 个正确端点，因此普通指标有 9 个可估数据，过程/贡献只有主队列 6 个可估数据。

| 指标 | 可估数据 | ELBO--指标 Spearman 中位数 | 方向一致 | winner 优于或等于正确端点中位数 |
|---|---:|---:|---:|---:|
| 观测点 NRMSE | 9 | -0.90 | 8/9 | 9/10 |
| dense NRMSE | 9 | -0.25 | 7/9 | 6/10 |
| loading error | 9 | -0.50 | 6/9 | 6/10 |
| feature total ISE | 9 | -0.20 | 6/9 | 8/10 |
| kernel relative ISE | 9 | -0.20 | 6/9 | 6/10 |
| covariance error | 9 | -0.20 | 6/9 | 6/10 |
| matched function recall | 9 | +0.30 | 6/9 | 8/10 |
| trajectory correlation | 9 | +0.60 | 7/9 | 8/10 |
| score correlation | 9 | +0.10 | 6/9 | 6/10 |
| factor-process NRMSE | 6 | -0.85 | 4/6 | 4/6 |
| complete-contribution NRMSE | 6 | -0.70 | 4/6 | 4/6 |

解释如下：

- 最大普通 ELBO 很可靠地偏向更好的观测点拟合，并在多数数据中同时偏向更好的 trajectory、process 和 contribution；
- 它对 dense、loading、feature、kernel 和 covariance 只有中等或不稳定的一致性；
- 某些数据的 winner 可以是某一科学指标上最差的正确端点，因此 ELBO 不是后验真值质量的万能排序；
- 但正式选择不允许使用真值指标。在当前候选中，ELBO 是合法、可复现且已证明能排除最严重错误/尺度端点的选择规则，不能因为它不是 oracle 就改用 NRMSE、ISE、R/P/L 或身份标签选解。

## 8. strict practical convergence 的含义

### 8.1 结构层面

| 端点组 | n | 200 sweep | strict converged | auxiliary endpoint pass | ever stable 5 | endpoint stable 5 |
|---|---:|---:|---:|---:|---:|---:|
| 全部 | 104 | 104 | 0 | 29 | 29 | 28 |
| 正确身份 | 45 | 45 | 0 | 8 | 7 | 7 |
| 非正确身份 | 59 | 59 | 0 | 21 | 22 | 21 |
| 冻结 winner | 10 | 10 | 0 | 3 | 3 | 3 |

严格收敛没有任何阳性端点，因而无法估计“strict 与科学质量”的直接关联。辅助 endpoint pass 的正确身份率为 8/29（27.6%），低于 fail 组的 37/75（49.3%）；endpoint-stable-5 组也只有 7/28 正确。辅助稳定性显然不能承担共享--特异身份判断。

### 8.2 在已经正确身份的端点内

条件在 45 条正确端点内，8 条 auxiliary endpoint pass 端点的科学指标中位数优于 37 条 fail 端点：

| 指标 | auxiliary pass | auxiliary fail |
|---|---:|---:|
| 观测点 NRMSE | 0.159 | 0.173 |
| dense NRMSE | 0.207 | 0.340 |
| loading error | 0.125 | 0.138 |
| feature total ISE | 0.071 | 0.138 |
| kernel relative ISE | 0.119 | 0.214 |
| covariance error | 0.420 | 0.667 |
| trajectory correlation | 0.980 | 0.970 |
| score correlation | 0.957 | 0.933 |
| factor-process NRMSE（主队列） | 0.231 | 0.271 |
| complete-contribution NRMSE（主队列） | 0.191 | 0.225 |

这说明末期数值稳定性在“身份已经正确”的条件下可能携带一定质量信息，但不能反向解释为成功判据：它漏掉 37/45 条正确端点，并接纳 21 条身份不正确端点；样本也不是随机分组。正确做法是继续把 strict convergence、objective eligibility 和科学输出稳定性分开报告，而不是把 strict gate 放宽或把 auxiliary pass 用作 winner 选择器。

## 9. 科学结论与决策

### 9.1 可以支持的结论

1. G12 在原始 6--9 时点场景的正确身份盆地是可达的，10 个数据均至少有一个正确端点；但单起点命中率波动大，多起点仍是必要预算。
2. 最大有效普通 ELBO 在这 10 个开发数据上提供了强的盆地支持：10/10 winner 正确、所有数据的最佳正确 ELBO 均超过最佳错误 ELBO。
3. 正确身份盆地内部大多数端点具有良好方向、feature/kernel coverage 和较好的过程相关，但存在少量尺度异常子盆地。
4. 当前多起点 ELBO 选择实质上管理了最危险的尺度异常：灾难性 dense/loading/covariance 尾部没有进入冻结 winner。
5. strict practical convergence 只描述末期增量，不描述结构或真值正确性；本批 0/104 strict 的结果不能用来否定科学上可用的 winner，也不支持仅为获得 strict 标记而全面延长。

### 9.2 不能支持的结论

- 不能把 45/104 当成未来数据的独立成功率，或把 10/10 winner 当成论文级 Monte Carlo 成功率；
- 不能宣称 G12 已彻底解决优化路径不稳定性；
- 不能宣称 weak separation 比 baseline 更容易；
- 不能把正确因子数、PPI 饱和或 auxiliary convergence 当成正确身份；
- 不能声称最大 ELBO 同时最优化 dense、loading、feature、kernel、covariance、process 和 contribution；
- 不能由本审计修改模型、先验、ELBO、初始化、退火、12 起点或默认 sweep 预算。

## 10. 下一步建议

1. 固化 Stage 5A v1 证据，继续把 G12 作为下一版开发协议的首选 R 初始化；保留历史初始化接口用于复现。
2. 保持 12 起点、一次预得分、Jaoua 退火和最大有效普通 ELBO 的 truth-free 选择；单起点命中率尚不足以删减多起点。
3. 不把 strict convergence 当作优化目标，也不因 0/104 strict 自动全面延长。若将来需要报告数值停止语义，继续单列 objective eligibility、terminal stability 和 scientific stability。
4. 若继续研究正确盆地可达性，下一阶段应是小规模、预注册的早期轨迹分流审计：在少量既有原始密度数据上，对正确、数量对但身份错、缺失因子三类起点追踪退火期与最初普通 sweep 的 block PPI、loading/score/process energy、RSS/ELBO 和共享--特异重叠。它先定位分流机制，不立即提出多种新算法。
5. 只有当上述轨迹证据反复显示同一关闭机制时，才比较一个可解释且等预算的优化日程候选；不得使用真值选择 winner，也不得同时混入多种技巧。
6. 与正确盆地可达性分开，后续论文级工作仍需独立验证维数选择、不同 p/信号强度、协变量 `d>0` 和更广泛场景；Stage 5A 不回答这些问题。

## 11. 关键机器结果

- `analysis_20260829_v1/AUDIT_QC.csv`
- `analysis_20260829_v1/AUDIT_COMPLETE.txt`
- `analysis_20260829_v1/INPUT_BINDING.csv`
- `analysis_20260829_v1/ALL_104_STANDARD_DENSITY_G12_ENDPOINTS.csv`
- `analysis_20260829_v1/EXACT_45_ENDPOINT_QUALITY.csv`
- `analysis_20260829_v1/FROZEN_10_WINNER_QUALITY.csv`
- `analysis_20260829_v1/PER_DATA_BASIN_SUPPORT.csv`
- `analysis_20260829_v1/CORRECT_BASIN_QUALITY_SUMMARY.csv`
- `analysis_20260829_v1/PER_DATA_EXACT_BASIN_METRIC_AUDIT.csv`
- `analysis_20260829_v1/WITHIN_DATA_ELBO_SCIENCE_ASSOCIATION.csv`
- `analysis_20260829_v1/ELBO_SCIENCE_ASSOCIATION_SUMMARY.csv`
- `analysis_20260829_v1/WINNER_WITHIN_EXACT_RANK_SUMMARY.csv`
- `analysis_20260829_v1/PRIMARY_PROCESS_CONTRIBUTION_BY_ROLE.csv`
- `analysis_20260829_v1/PRIMARY_PROCESS_CONTRIBUTION_SUMMARY.csv`
- `analysis_20260829_v1/PRIMARY_LOADING_SCALE_BY_BLOCK.csv`
- `analysis_20260829_v1/PRIMARY_LOADING_SCALE_SUMMARY.csv`
- `analysis_20260829_v1/CONVERGENCE_STRUCTURE_SUMMARY.csv`
- `analysis_20260829_v1/CONVERGENCE_CONDITIONAL_QUALITY_SUMMARY.csv`
