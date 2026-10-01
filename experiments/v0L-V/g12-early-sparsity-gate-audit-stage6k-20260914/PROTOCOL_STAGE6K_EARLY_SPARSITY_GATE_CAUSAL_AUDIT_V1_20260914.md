# Stage 6K：早期稀疏反馈/稠密门控配对因果审计

版本：V1，2026-09-14  
状态：开发性、已授权的小型机制审计；不是正式 Monte Carlo，也不是默认算法修改。

## 1. 问题与证据边界

在当前 Jaoua/Splinets 真值生成机制的困难数据 `s6g_minimal_L111_M2`
（数据种子 `84611001`）上，Stage 6H 的 12 个 G12 起点为 0/12 core；
Stage 6J 的 4 个配对 G12 起点也为 0/4 core。Stage 6I 的 oracle
局部化对照为 6/6 core，且正确端点的普通 T=1 ELBO 高于所有已观察到的
错误端点。因此该数据上存在更好的正确盆地，尚不能把失败归因于普通 ELBO
偏爱错误端点。

此前 2026-07-27 的小维度审计已经测试过 0/2/5/10 sweep 稠密门控；它只在
少量纯特异路径中有帮助，在混合模型中并非通用修复。因此本轮不是重新搜索
门控长度，而是在当前不同条件（G12、一次预得分、Jaoua 退火、`p=500`、
新真值生成器）下做一次有停止规则的因果复核。

## 2. 唯一干预与配对设计

- 数据：只复用 `s6g_minimal_L111_M2` 的 observation bundle，不生成新数据。
- 对照：只读继承 Stage 6J 的 4 个 G12 端点。
- 干预：新增 4 个 G12 端点，拟合种子与对照逐一配对：
  `84611101`--`84611104`。
- 拟合维度：`L_f=1`、`L_s=(1,1)`、所有 `M=2`、`K=5`。
- 共同路径：random + Gram 单位函数能量（G12）+ 1 次预得分 +
  Jaoua `c(1,1.9,100)` 退火 + 固定 400 个普通 T=1 sweep。
- 唯一改变：总 sweep 1--10 将共享和特异 loading inclusion coordinates
  暂时固定为 `gamma=1`；第 11 sweep 恢复原 spike-and-slab/PPI 更新。
- 每条 fit 内 `n_cpus=1`；4 条 fit 可在任务层并发。

`gate=10` 的理由仅为覆盖 Stage 6J 已观察到的早期塌缩窗口，并与历史审计
中最长的有界门控一致。不得根据本轮结果追加 2、5、20 等长度。

该干预不仅“推迟一个 omega 更新”：初始 gamma 也由普通路径的 0.5 暂时改为
1，并且 omega 会接收门控后的 gamma。因此准确名称是“短暂稠密门控/延迟
普通稀疏反馈生效”，不是先验或模型的永久修改。

## 3. 真值隔离、终点与选择

- 拟合、停止和真值盲选择只读取 observation bundle，不读取 sealed truth。
- 4 个新端点必须先全部终止并保存；主路径错误原样保留，不选择性重跑。
- objective eligibility 与 practical convergence 分开记录。
- 新干预臂内部仍按最大有效普通 T=1 ELBO 冻结一个 winner。
- 真值只在上述端点与 winner 冻结后用于科学评价。
- 不做跨维度选择、不做 continuation、不生成新 fit seed。

## 4. 记录指标

真值盲记录：

- 普通 T=1 ELBO、warning/error、运行时间、practical convergence；
- sweep 1、2、5、10、11、12、25、50、75、99、100、120、199、299、399、499
  的时间函数、score、loading driver、slab norm、unrestricted PPI 与实际 PPI；
- 门控是否在第 11 sweep 释放、释放后 sweep 数；
- PPI 因子数与观测点重构。

解封后记录：

- 载荷优先的一对一共享--特异身份（core）；
- observed/dense signal NRMSE、loading、feature、factor process、完整贡献和
  covariance error；
- 与同 seed Stage 6J G12 对照的配对差异。

## 5. 预注册判据与停止规则

主要判据是 4 个干预端点中的 core 命中数；Stage 6J 配对对照为 0/4。

- `>=2/4 core`，且正确端点普通 ELBO 不低于相应配对对照：支持“早期稀疏
  反馈是该困难数据的实质障碍”。这只允许把机制列为后续候选，不足以改默认。
- `1/4 core`：证据弱且不稳定，不改默认算法。
- `0/4 core`：关闭短暂稠密门控路线；接受该种子为 G12 的已知困难例，转回
  论文所需维度准则与含协变量模拟，不再追加门控长度或稀疏超参数搜索。

无论结果如何，本轮后停止并报告。Stage 6K 不修改模型、先验、CAVI、ELBO、
公开接口、G12 默认路径或任何历史冻结实验。
