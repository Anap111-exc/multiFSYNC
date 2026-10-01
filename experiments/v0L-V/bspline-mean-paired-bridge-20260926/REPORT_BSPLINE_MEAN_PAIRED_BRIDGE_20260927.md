# B 样条均值小型配对桥接：科学评价报告

实验 ID：`BSPLINE_MEAN_PAIRED_BRIDGE_V1_20260926`  
评价日期：2026-09-27  
性质：开发阶段小型桥接，不是论文级 Monte Carlo。

## 1. 执行与真值隔离

- 2 个新逻辑数据种子、2 个配对均值臂、每臂每种子 12 个配对起点，共 48 条 fit。
- 48/48 完成，0 error，48/48 objective eligible。
- 每个“数据种子 × 均值臂”内的最大有效普通 ELBO winner 在解封前已经冻结；禁止跨均值臂比较 ELBO。
- 用户于 2026-09-27 授权解封后，48/48 endpoint 和 4/4 winner 完成科学评价，0 评价错误、0 评价 warning。
- 未启动新 fit、continuation、racing 或 800-sweep 延长。

## 2. 正确盆地可达性与 ELBO 选择

| 数据 | 均值臂 | core 命中 | 因子数 1/1/1 | exact structure | strict practical | winner |
|---|---|---:|---:|---:|---:|---|
| meanb_01 | 解析周期 P | 1/12 | 2/12 | 1/12 | 10/12 | start 08，core |
| meanb_01 | B 样条 B | 1/12 | 3/12 | 1/12 | 11/12 | start 08，core |
| meanb_02 | 解析周期 P | 4/12 | 8/12 | 4/12 | 11/12 | start 07，core |
| meanb_02 | B 样条 B | 6/12 | 8/12 | 6/12 | 11/12 | start 07，core |

合计：P 为 5/24 core，B 为 7/24 core。24 组相同拟合种子的配对转移为：

- 两臂均 core：5；
- P core、B 非 core：0；
- P 非 core、B core：2；
- 两臂均非 core：17。

因此在本次实现中没有发现 B 样条均值系统性降低可达性的证据；它净增加 2 个 core endpoint。但数据种子只有 2 个，不能把 5/24 与 7/24 的差异解释为稳定性能提升。

四个 ELBO winner 均为 core 且为 exact structure。meanb_01 每臂只有一个 core，winner 相对第二名的 ELBO 间隔分别为 8814.68（P）和 8807.06（B）；meanb_02 存在多个 core，winner 相对第二名仅差 14.30（P）和 26.08（B），但最佳 core 相对最佳非 core 仍分别高 7486.21 和 7480.10。这表明当前两种均值真值下，普通 ELBO 都能区分正确盆地和错误盆地；meanb_01 的 1/12 命中同时说明优化可达性仍然脆弱。

## 3. 均值函数恢复

四个冻结 winner 的结果为：

| 数据 | 均值臂 | observed NRMSE | dense NRMSE | dense ISE | 投影下界 ISE | excess ISE |
|---|---|---:|---:|---:|---:|---:|
| meanb_01 | P | 0.13468 | 0.13916 | 0.004357 | 0.00000623 | 0.004351 |
| meanb_01 | B | 0.12805 | 0.13320 | 0.003992 | 0.00001528 | 0.003976 |
| meanb_02 | P | 0.13750 | 0.13997 | 0.004408 | 0.00001293 | 0.004395 |
| meanb_02 | B | 0.13436 | 0.13614 | 0.004170 | 0.00001025 | 0.004160 |

B 臂的 dense NRMSE 在两个种子中均略低于配对 P 臂，但差值只有 -0.0060 和 -0.0038。所有均值投影下界仅约 `6.2e-6` 至 `1.53e-5`，不足总均值 ISE 的 0.4%；误差主要是有限样本估计误差，而不是 B 样条真值落在拟合空间之外。observed 与 dense NRMSE 很接近，也没有点间或边界爆炸。

## 4. 冻结 winner 的完整科学质量

按两个数据种子平均：

| 指标 | P | B |
|---|---:|---:|
| observed signal NRMSE | 0.13266 | 0.13243 |
| dense signal NRMSE | 0.21734 | 0.21687 |
| feature total ISE | 0.06432 | 0.06329 |
| kernel relative ISE | 0.09759 | 0.09695 |
| covariance operator relative error | 0.39340 | 0.39311 |
| loading relative L2 error | 0.11338 | 0.11351 |
| loading absolute cosine | 0.999758 | 0.999759 |
| FPCA score absolute correlation | 0.96887 | 0.96979 |
| factor trajectory absolute correlation | 0.98812 | 0.98877 |
| factor-process NRMSE | 0.18166 | 0.17722 |
| complete-contribution NRMSE | 0.15276 | 0.14878 |

所有 winner 均满足：

- 正确共享/研究 1 特异/研究 2 特异身份，无 missing、misplaced、duplicate 或 extra；
- 选择因子数为 1/1/1；
- 三个因子的有效 M 均为 2；
- 8/8 真 FPCA 成分均被覆盖。

科学质量的主要差异来自数据种子，而不是均值臂。例如 dense signal NRMSE 在 meanb_01 约 0.166--0.167，在 meanb_02 约 0.267；协方差算子误差相应约为 0.272--0.274 和 0.513--0.514。由于每个逻辑种子内因子、载荷、得分、噪声和观测设计均完全配对，P/B 几乎相同的结果说明更换均值生成方式没有改变残差共享—特异结构的主要科学结论。

在 5 对“两臂均 core”的同起点 endpoint 中，B-P 的平均差异为：均值 dense NRMSE -0.00621、observed signal NRMSE -0.000175、dense signal NRMSE -0.000338、载荷相对 L2 -0.00117、过程 NRMSE -0.00322。其余指标差异也很小；这些数字只作配对描述，不作显著性或优势声明。

## 5. 收敛、warning 与计算成本

- strict practical convergence：B 22/24，P 21/24，总计 43/48。
- 5 条非严格收敛轨迹均保留且仍为 objective eligible；它们均不是 core。
- 但 43 条严格收敛轨迹中有 31 条仍不是 core，因此严格收敛远不足以识别正确盆地。
- 每条 fit 平均约 36.0--36.2 分钟，范围约 31.5--38.6 分钟；两臂计算成本基本相同。
- 8 路并发从正式启动至 supervisor 完成约 3 小时 41 分。预检的 2.24 小时低估了实际成本。

## 6. 结论与边界

本桥接支持以下结论：

1. `independent_splinet` B 样条均值与当前 multiFSYNC 拟合流程兼容；没有发现实现、投影或边界异常。
2. B 样条均值没有恶化正确盆地可达性；四个真值盲 winner 全部正确。
3. 残差共享—特异结构、FPCA 成分、载荷、得分和过程恢复在两臂间近乎不变。
4. 因此可以在下一轮含协变量开发实验中把独立 B 样条均值设为主场景，并保留解析周期均值作为配对敏感性场景；不需要重跑先前所有历史开发实验。
5. 本实验只有 2 个数据种子、强识别载荷、`d=0`、真实 L/M 和周期 Splinets 均值，不能宣称 B 样条均值普遍更优，也不能消除优化可达性问题。
6. 当前 B 样条均值仍采用周期边界，以公平隔离“解析式 vs 随机 B 样条生成”。它不等于已经验证非周期 B 样条均值；若论文时间域明确无周期含义，应把该边界设定写清楚，而不是把本桥接外推为非周期验证。

## 7. 机器结果位置

- 完整评价归档：`bspline_mean_bridge_evaluation_20260927_v1.tar.gz`
- 解压机器结果：`bspline_mean_paired_bridge_20260926_v1_evaluation_20260927/`
- 关键表：
  - `STRATUM_BASIN_AND_WINNER_SUMMARY.csv`
  - `WINNER_SCIENTIFIC_PRESENTATION.csv`
  - `MEAN_ARM_SCIENCE_SUMMARY.csv`
  - `PAIRED_CORE_TRANSITION_SUMMARY.csv`
  - `PAIRED_START_ENDPOINT_TRANSITIONS.csv`
  - `MEAN_PROJECTION_FLOOR_BY_DATA_ARM.csv`

归档 SHA256：`f42e58f83a4dbf85fb53973bfd2da630428e3977e010e3b89d590735601800a8`。
