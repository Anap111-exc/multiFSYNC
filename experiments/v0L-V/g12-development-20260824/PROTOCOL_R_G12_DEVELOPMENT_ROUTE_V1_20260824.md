# R-G12 下一版开发路线协议 v1

日期：2026-08-24  
协议标识：`R_G12_DEVELOPMENT_ROUTE_V1_20260824`  
状态：`approved_for_next_development_protocol`  
正式实验执行授权：`FALSE`

## 1. 决策

基于两轮独立的新数据验证，将 `gram_unit_energy` 提升为下一版 post-v0L-V 开发协议的首选 R 初始化。该决策只改变协议显式传入的随机函数初始化几何，不修改统计模型、先验、CAVI、完整 ELBO、预得分更新、退火计划、主预算或评价公式。

公共接口 `bayesSYNC_multi_pre_score()` 的历史默认值继续为 `current_1_over_m`。新协议必须显式传入：

```r
function_initialization = "gram_unit_energy"
```

不得通过静默更改包默认参数，使历史脚本在未声明的情况下改用 G12。

## 2. 证据基础

### 2.1 G12 精简独立确认

- 6 个预注册新数据；
- 每方法每数据 12 个配对起点；
- G12 严格正确起点 27/72，current 为 8/72；
- G12 正确 truth-free winner 6/6，current 为 5/6。

### 2.2 G12 跨场景推广确认

- 6 个预注册新数据，覆盖基线强可辨识、弱载荷分离和稀疏观测；
- 每方法每数据 8 个配对起点；
- G12 严格正确起点 23/48，current 为 2/48；
- G12 正确 truth-free winner 6/6，current 为 2/6。

两轮描述性合并：G12 为 50/120 严格正确起点和 12/12 正确 winner；current 为 10/120 和 7/12。该合并不是论文级 Monte Carlo，但已达到把 G12 提升为下一版开发路线的证据要求。

## 3. 下一版 R 路线

路线标识：

`random__gram_unit_energy__pre1__jaoua_default__multistart12`

固定配置：

- `initialization="random"`；
- `function_initialization="gram_unit_energy"`；
- Gram/L2 单位函数能量校准只作用于迭代零时刻的随机函数变分均值；
- 1 次有界预得分；
- Jaoua 默认退火 `c(1,1.9,100)`，即 99 个退火 sweep；
- 退火后普通 `T=1` CAVI；
- 主验证路线使用 12 个预注册起点；
- 旧 strict practical convergence 若连续满足则允许提前停止，否则最多 200 个普通 T=1 sweep；
- objective eligibility 与 strict practical convergence 分开；
- 在 objective eligible 且普通 T=1 ELBO 有限的终点中，按最大普通 T=1 ELBO 选择唯一 winner；
- tie 只按预注册的 `fit_id` 字典序处理；
- 每 fit `n_cpus=1`，并行只发生在独立 fit 层；
-不进行实际 racing；
- 不自动 continuation；
- 不使用真值、结构标签、R/P/L、NRMSE 或 ISE 进行停止或选择。

小型开发诊断若因预算使用少于 12 个起点，必须在独立协议中预注册；这不改变 R-G12 主路线的 12 起点定义。

## 4. 历史复现路线

`current_1_over_m` 保留为历史复现和消融对照选项：

```r
function_initialization = "current_1_over_m"
```

它不再是下一版开发协议的首选初始化。不得把 current 的起点和 G12 的起点合并成一个 24 起点正式候选池，也不得在看到真值后选择使用哪条初始化路线。

下列历史文件和结果保持原样，不做追溯修改：

- R-only v1/v2/v3/candidate2 协议；
- v0L-V 总协议 v1/v2/candidate2；
- candidate2 manifest 和冻结绑定；
- 已完成的 v0J、v0K、v0K+、v0L、v0L-V 及两轮 G12 实验。

## 5. 暂不改变的部分

- multiFSYNC 模型定义；
- 先验及超参数；
- CAVI 更新公式；
- 完整 ELBO；
- 预得分的更新顺序和次数；
- Jaoua 退火；
- 200-sweep 主预算；
- 最大有效普通 ELBO winner 规则；
- 载荷优先、过程/贡献二阶段确认的评价框架；
- `lambda_orth=0` 默认值。

## 6. convergence 边界

旧 strict practical convergence 暂时只保留为计算停止记录，不是 objective eligibility、winner 排序或科学成功门槛。2026-08-24 的离线语义审计独立进行；在该审计完成前不修改阈值，也不全面运行 200→400/800 continuation。

若后续提出新的停止语义，必须：

1. truth-free；
2. 与 objective eligibility 分开；
3. 先对既有轨迹离线 replay；
4. 只在无法由既有轨迹判断时，另行预注册极小 continuation 审计；
5. 不利用本批真值逐阈值调参后在同批数据上宣称验证成功。

## 7. 执行授权边界

本协议确认路线方向，但不授权任何新实验。建立具体数据/fit manifest、生成数据、启动拟合或 continuation，仍需用户对相应版本化实验协议单独授权。

