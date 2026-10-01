# Stage 6J：真值盲可达性转折点与 residual-FPCA 配对审计 v1

协议 ID：`STAGE6J_TRUTH_FREE_REACHABILITY_TRANSITION_V1_20260914`  
性质：开发期、单数据种子、机制定位；不是正式论文 Monte Carlo。  
数据：只读复用 Stage 6G 的 `s6g_minimal_L111_M2` 观测包；拟合阶段不得读取 sealed truth。

## 1. 问题与边界

Stage 6I 已在同一数据、真实固定维度 `L_f=1, L_s=(1,1), M=2` 下证明：正确结构盆地存在，且其普通 `T=1` ELBO 高于已产生的错误随机端点；但现行随机 Gram 单位函数能量路径没有生成正确端点。Stage 6J 只回答：

1. 特异块在预得分、退火还是最初普通 `T=1` sweep 中关闭；
2. 一个真值盲、study-aware 的 residual-FPCA 初始化能否在相同种子和预算下提高候选盆地可达性；
3. 若仍失败，下一步应优先检查初始化、预得分/退火时序，还是早期稀疏反馈。

本阶段不修改模型、先验、普通 ELBO、CAVI、因子维度、预得分次数、退火温度或停止规则；不测试跨维度选择准则，也不运行 800 sweep。

## 2. 两个配对拟合臂

共同设置：

- 数据种子：`84611001`（既有观测包）；
- 拟合种子：`84611101:84611104`，两臂严格配对；
- `L_f=1, L_s=(1,1), M_f=2, M_s=((2),(2)), K=5`；
- 一次预得分；Jaoua 退火 `c(1,1.9,100)`；
- 99 个退火 sweep + 400 个普通 `T=1` sweep；
- `lambda_orth=0`，每 fit `n_cpus=1`；
- 外层最多 8 个独立 fit 并发；
- endpoint eligibility 与 practical convergence 分开记录。

拟合臂：

1. `random_gram_unit_energy`：`initialization="random"`，随机函数变分均值经 Gram 单位能量校准。这就是历史简称 G12 的拟合初始化。
2. `residual_fpca`：`initialization="residual_fpca"`，仅用观测数据构造 study-aware 谱初始化；保持同一次预得分、同退火和同 CAVI。`function_initialization` 字段对该臂不参与随机函数校准。

若 residual-FPCA 某块因数值秩不足触发代码既有的 seeded-random fallback，fit 仍原样保留并记录 warning、失败块和 fallback 次数；它不被伪装成完整 residual-FPCA，也不选择性重跑。

总计 8 条新 fit。Stage 6H/6I 的端点只作背景证据，不替代本次配对拟合。

## 3. 转折点留痕

`trace_sweeps` 是从 1 开始的总 CAVI sweep 编号，而不是普通 `T=1` sweep 编号。预注册记录：

- 总 sweep `1:120`：覆盖全部 99 个退火 sweep及最初 21 个普通 `T=1` sweep；
- 总 sweep `129,149,199,299,399,499`：对应约普通 `T=1` sweep `30,50,100,200,300,400`。

每个记录点保留共享/特异载荷、时间函数和得分的 driver、precision、PPI、posterior-mean norm 等既有诊断字段。诊断只读取更新中已有量，不改变更新。

## 4. 真值隔离与终点选择

- prepare、smoke、fit、停止、eligibility、汇总和 winner 选择均不得读取 truth bundle；
- 同一固定维度的 8 个有效端点按最大普通 `T=1` ELBO作唯一 truth-free winner；
- warnings 不得被 `suppressWarnings()` 丢弃；error 原样保留，不选择性重跑；
- 真值科学评价需用户另行明确授权，并写入新评价目录。

## 5. 预注册判断与后续动作

解封后的主要判断单位是 4 对拟合种子，而不是单条最好结果：

1. residual-FPCA 若产生至少 `2/4` core，且其正确端点普通 ELBO不低于错误盆地：把它列为候选真值盲初始化，再用 2 个新数据种子确认；不立即替换开发默认路径。
2. residual-FPCA 初始结构证据存在、但在退火中与 G12 同期关闭：优先做一次小型早期稀疏反馈/退火时序审计，不再堆叠初始化补丁。
3. residual-FPCA 在预得分前后即丢失 study-specific 证据：检查 shared-first residual 分解和预得分更新顺序。
4. 两臂都能产生正确端点，而最大普通 ELBO仍选错：才把问题升级为固定维度端点选择/目标问题。
5. 两臂均失败：停止宣称初始化问题已解决；下一阶段只允许一个预注册的 `Beta(1,p)` 或稀疏门控时序敏感性实验。

Stage 6J 只有一个已开发数据种子，任何正结果都只是机制证据，不能宣称跨数据推广。
