# Stage 7C：Route-B 含协变量独立确认协议 v1 草案

协议 ID：`COVARIATE_ROUTEB_STAGE7C_V1_20260927`  
状态：`draft_not_authorized`  
日期：2026-09-27  
计划规模：3 个新逻辑数据种子 × 3 个配对臂 × 12 个拟合起点 = 108 fits

## 1. 目的与证据定位

Stage 7C 独立确认 36-fit 开发先导的三个结论：

1. 零真实效应时，含协变量模型不会系统地产生实质性伪 beta；
2. 非零跨研究共享 B 样条 beta 能被恢复，同时残差共享—特异因子结构仍可识别；
3. 对同一响应遗漏协变量时，beta 信号会污染均值或因子块，并损害残差结构恢复。

本阶段不比较初始化方法，不选择 `K`，不进行 `L_f/L_s` 或 `M_f/M_s` 维度选择，不修改模型、先验、CAVI、ELBO、停止规则或评价公式。它是三个独立数据重复的确认实验，仍不是论文最终 Monte Carlo。

36-fit 先导的逻辑数据种子 `covrbp_01` 不进入 Stage 7C；其结果只用于冻结设计和给出开发性警戒阈值，不参与 Stage 7C 成功比例的分母。

## 2. 三臂配对设计

每个新逻辑数据种子 `covrb7c_01`、`covrb7c_02`、`covrb7c_03` 均生成以下三臂：

| 臂 | 数据生成 | 拟合输入 | 配对关系 |
|---|---|---|---|
| C0 | `d=2`，保存 Z，真实 `beta=0` | Y、Z、time | 与 C1 共用潜在真值、时点、Z 和噪声，仅移除 beta 贡献 |
| C1 | `d=2`，非零跨研究共享 B 样条 beta | Y、Z、time | 正确调整臂 |
| C2 | 直接复用同一逻辑种子的 C1 响应和时点 | Y、time，`d_fit=0` | 与 C1 的观测数据逐元素相同，仅遗漏 Z |

C0/C1 必须逐元素核验：除 beta 真值、完整 signal 和响应外，其均值、特征函数、得分、载荷、时间、Z 与噪声完全相同。C1/C2 必须逐元素核验 Y 和 time 完全相同。C2 不得重新生成响应或噪声。

同一逻辑数据种子和同一起点编号在 C0/C1/C2 使用相同 `fit_seed`。因为 C2 的参数维数不同，这表示预注册的配对起点标签，不声称随机数流逐项相同。

## 3. 固定数据生成与拟合配置

- `S=2`，`n_s=(30,30)`，`p=500`；
- 生成协变量维数 `d=2`，`Z_s[i,r] iid N(0,1)`，不与因子得分人为正交化；
- 真实 `L_f=1`、`L_s=(1,1)`；
- `M_f=(2)`、`M_s=((2),(2))`，FPCA 得分方差权重为 `1/m^2`；
- 每人 6--9 个独立 `Uniform(0,1)` 不规则观测时点；
- `sigma_eps=0.3`，`mean_amp=0.6`，`n_dense=201`；
- 特征函数：`iid_common_splinet_stiefel`，公共维数 8、非周期 regular Splinets、零边界、平滑权重 3；
- 均值：`independent_splinet`，周期三次 Splinets、公共维数 8、平滑参数 `1.519826404906808`；
- 载荷：`identified_loadings=TRUE`、`bool_sparse_loadings=TRUE`、`prop_sparse=0.9` 的受控强识别场景；
- 生成器缓存和拟合均使用 `K=9`；`lambda_orth=0`；
- C1 beta：每个 `(j,r)` 使用独立三次普通 B 样条，内部结点 `0.25,0.5,0.75`、边界 `0,1`、`intercept=TRUE`、基维数 7；各曲线按 201 点梯形积分匹配参考 `0.3*sin(pi*r*t+phase_jr)` 的连续均值与中心化 `L2` RMS；在真实观测时点直接求值，不由 dense grid 线性插值；
- C0 beta 恒为零；C2 复用 C1 的非零 beta 数据但拟合不接收 Z。

本阶段保持当前非周期零边界特征函数、周期 B 样条均值和受控强识别载荷，不混入 Jaoua 式随机载荷、自由边界特征函数或正弦 beta，以避免同时改变多个生成轴。

## 4. 待用户确认的预注册 seeds

| 逻辑数据 | 基础数据/时点/得分/Z/噪声 | 特征函数 | 载荷 | 均值 | beta | 12 个配对拟合 seeds |
|---|---:|---:|---:|---:|---:|---:|
| `covrb7c_01` | 99711001 | 99712001 | 99713001 | 99714001 | 99715001 | 99716101--99716112 |
| `covrb7c_02` | 99711002 | 99712002 | 99713002 | 99714002 | 99715002 | 99716201--99716212 |
| `covrb7c_03` | 99711003 | 99712003 | 99713003 | 99714003 | 99715003 | 99716301--99716312 |

这些整数已在当前仓库的非结果源文件、协议和 manifests 中做精确冲突检索，未发现历史复用。它们在用户明确确认前仅是候选，不得用于生成数据。

## 5. 固定拟合路径与计算资源

每条 fit 固定：

- random 初始化；
- `gram_unit_energy`（G12）函数变分均值初始化；
- 1 次有界预得分；
- Jaoua 退火 `c(1,1.9,100)`，即 99 个退火 sweep；
- 随后固定 400 个普通 `T=1` sweep，`maxit=499`；
- `K=9`、`n_g=51`、`lambda_orth=0`；
- 拟合维度与真值相同：`L_f=1`、`L_s=(1,1)`、所有 `M=2`；
- 每条 fit 内部 `n_cpus=1`；
- ECS 外层最多并发 8 条相互独立 fit；
- 不使用 racing、自动 continuation、400→800 延长或真值辅助初始化。

每条 fit 独立目录至少保存完整 fit RDS、普通 ELBO 与分量、warning/error、运行时间、诊断轨迹以及 `COMPLETE/ERROR` 标记。error 或 objective-ineligible 记录原样保留，不选择性重跑。

## 6. Objective eligibility 与真值盲选择

objective eligibility 与 strict practical convergence 分开。候选需满足：

- 完整执行预注册退火和 400 个普通 `T=1` sweep；
- 最终温度为 1；
- 最终普通 ELBO 及其分量均为有限值；
- 普通目标有效；
- 普通阶段 ELBO decrease 计数为 0；
- 普通阶段 solver jitter 计数为 0。

在每个“逻辑数据种子 × 臂”内部，从 objective-eligible 候选中按最大最终普通 `T=1` ELBO 选一个 winner；并列时按 `fit_id` 字典序。共冻结 9 个 winner。禁止跨 C0/C1/C2 或跨数据种子比较 ELBO。

所有 108 条拟合的停止、eligibility 和 winner 选择均不得访问真实结构、beta、NRMSE、ISE、R/P/L 或其他真值字段。9 个 winner 和选择输入哈希冻结后停止，等待用户另行授权解封真值。

## 7. 解封后的固定科学评价

只有用户单独授权后才评价全部起点和 9 个 truth-free winner。评价沿用 36-fit 先导的既定实现：

- C0：伪 beta 绝对 ISE、总能量、相对对应 C1 真 beta 的能量比及伪贡献 RMS；零真值下不报告以零为分母的 relative ISE；
- C1：逐 beta ISE/relative ISE/连续 `L2` cosine/能量比，以及观测点和 dense-grid beta 完整贡献误差；
- C2：与同 seed、同 start 的 C1 配对比较，评价 beta 被均值或因子块吸收的程度；
- 三臂：包含 `mu + Z beta + factors` 的完整 observed/dense signal 指标、均值、载荷、特征函数投影下界与 excess ISE、得分、过程、完整贡献、协方差核、missing/misplaced/duplicate/extra、PPI、ELBO、运行时间、warning 和 strict convergence；
- 结构身份先由载荷方向与参数块一对一匹配，再评价连续函数与过程；
- strict convergence 只作优化状态描述，不作为科学正确性的替代标签。

旧的 `d=0` signal 指标因遗漏 `Z beta`，不得用作 C0/C1/C2 的完整 signal 主指标。

## 8. 预注册判据与停止规则

### 8.1 主要确认判据

三个独立 C1 数据池中：

1. 每个 12-start 候选池均至少出现 1 个 core/exact 端点；否则记录为 12 起点可达性失败，不追加起点；
2. 至少 2/3 个最大 ELBO winner 达到 core、exact 且 8/8 true feature coverage；
3. 至少 2/3 个 C1 winner 满足 beta relative ISE `<=0.20` 且平均连续 `L2` cosine `>=0.90`。

C0 中，至少 2/3 个 winner 的估计 beta 总能量不超过同 seed C1 真实 beta 总能量的 10%。

C2 作为遗漏协变量因果对照，不要求必然选择错误因子数；预注册方向是：至少 2/3 个 seed 的 C2 winner observed full-signal NRMSE 高于配对 C1 winner，且 C2 在共享—特异身份或完整贡献中表现更差。若不满足，应报告遗漏效应可被当前模型吸收的边界，而不是修改评价口径。

上述阈值是从单 seed 开发先导冻结出的确认门槛，不是普遍一致性或论文最终性能声明，也不得参与拟合选择。

### 8.2 失败与异常处理

- 数据配对 QC、密封、manifest 或源绑定失败：整批在拟合前停止；
- 单条 fit error/ineligible：保留原记录，不重跑；
- 某 seed 的 12 个候选全部 ineligible：该 seed-arm 无 winner，原样报告；
- 某 C1 池无 core：完成既定评价后只做只读归因，不自动增加起点、延长或调参；
- 任何协议、源码或 seed 变更均需新版本文件，不覆盖本草案。

## 9. 时间预算

36-fit 先导累计约 37.05 fit-hours，其中 C0/C1 每条约 76 分钟、C2 每条约 38 分钟。按相同配置线性外推，Stage 7C 约 111 fit-hours。8 核 ECS 每 fit 单核、最多 8 个外层 worker，考虑长尾、数据生成、密封、汇总和调度损耗，预留约 15--18 小时墙钟时间；科学评价通常为分钟级。

如基准测得单 sweep 明显偏离先导，不自动改变预算或 sweep 数，只报告并等待用户决定。

## 10. 当前授权边界

截至本草案创建时：

- `protocol_confirmed = FALSE`；
- `seed_schedule_confirmed = FALSE`；
- `execution_authorized = FALSE`；
- `stage7c_data_generated = FALSE`；
- `stage7c_fits_started = FALSE`；
- `truth_unseal_authorized = FALSE`；
- `scientific_evaluation_started = FALSE`；
- ECS 不因本草案而启动。

用户确认协议和 seed 表后，下一步才是实现最小 runner、运行本地解析/定向测试/微型 smoke；这些通过并报告后，才可启动 ECS 正式 Stage 7C 拟合。
