# Route B 非周期真值配对拟合桥接协议 v1

日期：2026-09-22  
状态：已获用户授权执行；开发性桥接，不属于正式 v0L-V 或论文级 Monte Carlo。

## 1. 目标

验证新的非周期 `iid_common_splinet_stiefel` 真值机制能否被 multiFSYNC
稳定恢复，并配对比较拟合空间 `K=5/8/9`：

- `K=5`：原低秩配置，作为表示能力不足对照；
- `K=8`：1000 次无拟合审计选出的最小通过值；
- `K=9`：投影下界更低的相邻敏感性值。

本实验不重新选择初始化、退火、CAVI、先验或 ELBO，也不跨 K 比较 ELBO。

## 2. 数据与真值隔离

- 2 个新逻辑数据种子；
- `S=2`, `n_s=(30,30)`, `p=500`, `d=0`；
- `L_f=1`, `L_s=(1,1)`，所有真 `M=2`；
- 每位受试者 6--9 个不规则时点；
- `sigma_eps=0.3`, `mean_amp=0.6`；
- 真值函数：非周期 `iid_common_splinet_stiefel`，`J_truth=8`,
  `smoothness_alpha=3`；
- 载荷：受控、强可辨识、固定 90% 稀疏，并按既有 v0F 规则校准；
- Jaoua 式随机载荷不与本桥接混合，留作独立敏感性实验。

生成器中的 `K=8` 只用于生成对象内部的 O'Sullivan 设计缓存；响应真值函数由
外部 Splinets 机制产生。每个逻辑数据种子只生成一份 `Y/Z/time_obs`，随后
`K=5/8/9` 拟合共享完全相同的观测数据。

拟合、停止与 winner 选择阶段不得读取密封真值。只有全部拟合终止且每个
`data_id x K` 层内的 winner 已冻结后，才可在用户另行授权下解封科学评价。

## 3. 配对起点与拟合路线

- 每个 `data_id x K` 使用 8 个配对随机起点；
- 相同 `start_index` 在两个数据种子和三个 K 下使用相同 `fit_seed`；
- 总拟合数：`2 x 3 x 8 = 48`；
- random 初始化；
- Gram 单位函数能量初始化（G12）；
- 1 次预得分；
- Jaoua 默认退火 `c(1,1.9,100)`，共 99 个退火 sweep；
- 普通 `T=1` CAVI 固定运行 400 sweep；
- 每条 fit `n_cpus=1`，BLAS/OpenMP 线程均固定为 1；
- ECS 最多 8 条独立 fit 并发；
- 不 racing、不 continuation、不自动 800 sweep、不选择性重跑失败条目。

增加起点数与并发数的理由是：8 个起点比原草案 6 个提供更高的盆地覆盖，且
与 8 核服务器一一对应；每条 fit 的统计路线和计算预算不变。

## 4. 真值盲终点选择

在每个 `data_id x K` 层内，从 objective-eligible、普通 `T=1` ELBO 有限的
终点中选最终普通 ELBO 最大者；并列时按 `fit_id` 字典序。禁止：

- 跨 K 比较普通 ELBO；
- 使用结构标签、ISE、NRMSE、载荷、R/P/L 或任何真值指标选 winner；
- 因单条失败而只重跑该条。

## 5. 启动门槛与资源控制

正式桥接批次前依次执行：

1. runner 解析和 manifest 自检；
2. 非正式微型 smoke，覆盖 K=5/8/9 的生成、拟合和真值盲选择；
3. 在全尺寸数据上做 K=9 的短 sweep 计时/内存预检；
4. 仅当没有错误、没有交换空间压力且峰值内存支持 8 并发时启动 48 条批次。

若预检显示 8 并发不安全，只降低独立 fit 并发数，不改变单条 fit、起点数或
科学设计。

## 6. 后续科学评价（本次启动阶段不执行）

另行授权解封后，评价所有端点及真值盲 winner：

1. 正确共享/特异核心、missing/misplaced/duplicate/extra；
2. 观测点与 dense-grid signal NRMSE；
3. 特征函数总 ISE、拟合空间投影下界和 excess ISE；
4. 协方差核/算子误差；
5. 载荷方向、FPCA 得分、因子过程和完整贡献；
6. ELBO 盆地、跨起点独立支持、运行时间、warning 和 practical convergence。

判定规则：若 `K=8/9` 在两个新种子上都不能产生正确核心端点，停止扩展并归因
可达性；若 K=8 与 K=9 科学恢复相近，优先 K=8；只有 K=9 在 dense、特征函数
或协方差恢复上出现稳定、实质改善且成本可接受时，才考虑 K=9。

## 7. 冻结边界

这是开发桥接，不修改历史实验，不启动正式 v0L-V，不修改模型、先验、CAVI、
ELBO、G12、退火或停止语义。执行结束前保持
`truth_unseal_authorized=FALSE` 和 `scientific_evaluation_started=FALSE`。
