# 既有轨迹 convergence 语义离线审计协议 v1

日期：2026-08-24  
协议标识：`G12_CONVERGENCE_SEMANTICS_OFFLINE_AUDIT_V1_20260824`  
状态：`authorized_read_only_audit`  
新拟合授权：`FALSE`  
continuation 授权：`FALSE`

## 1. 目的

只读合并两轮已经完成并解封评价的 G12 实验，检查旧 strict practical convergence 的各门是否具有合理的计算停止语义，以及 strict 标记是否能代表科学输出稳定或正确盆地。

本审计不寻找使真值恢复率最大的阈值，不实现新停止规则，也不修改已完成 fit。真值只在拟合和 truth-free winner 冻结后用于 post-hoc 语义检验。

## 2. 数据来源

审计使用两轮分析脚本从完整 `fit$practical_diagnostics` 轨迹重建的 endpoint/streak 表，并与相应全终点科学评价按 `fit_id` 精确连接：

- G12 精简独立确认：144 条终点；
- G12 跨场景推广确认：96 条终点；
- 合计 240 条，current 与 G12 各 120 条；
- 共 12 个预注册新数据。

复制进本审计目录的 source CSV 是只读快照；不重新读取或修改历史 RDS。

## 3. 固定旧 strict 规则

沿用 candidate2 的既有阈值，不在本审计中改变：

- `min_t1=80`；
- 短窗口 20、长窗口 60；
- 连续通过 5 次；
- objective rate gate；
- fitted NRMSE：短窗 `1e-3`、长窗 `3e-3`；
- RSS relative change：短窗 `1e-3`、长窗 `6e-3`；
- variable-PPI 最大绝对变化：短窗和长窗均 `1e-2`。

既有分析还记录两种纯离线替代语义：

- 用 variable-PPI 第 95 分位数和 factor-PPI 最大变化替代 elementwise PPI max；
- 在上述替代上进一步移除 objective gate，作为 output-only 诊断。

这些替代只用于审计，不自动成为新停止规则。

## 4. 预先确定的问题

1. strict convergence 在两方法、两轮实验中的发生率；
2. strict 与严格正确 `1/1/1` 的 post-hoc 交叉关系；
3. objective、fitted、RSS、PPI-max 在终点的阻断频率；
4. 去掉 PPI-max 是否足以连续通过 5 次；
5. quantile+factor-PPI 替代是否足以连续通过 5 次；
6. 各规则在正确与错误终点中的通过率，判断其是否明显偏好塌缩/低活动盆地；
7. strict 与重构、conditional feature、载荷、得分和 kernel 指标的关系；
8. 现有证据是否足以支持全面延长到 400/800。

## 5. 解释规则

- strict 是计算停止事件，不把它当作科学分类器；所谓 sensitivity/false-positive 只用于检查语义关联；
- 起点不是独立数据重复，不能把 240 条终点当作论文级样本；
- conditional feature ISE 必须与 matched/total coverage 同时解释；
- 若 PPI-max 普遍失败但 quantile/factor 替代仍无法覆盖科学正确终点，则不得宣称“只改 PPI 即解决 convergence”；
- 若正确终点系统性比错误/塌缩终点更难通过 fitted/RSS，需检查停止量是否受到活动维数、旋转或缓慢幅度更新影响；
- 本审计不能证明 400/800 没有收益，因为没有新增 continuation；它只能判断是否已有充分理由全面延长。

## 6. 输出和停止条件

输出版本化 CSV、QC、sessionInfo 和报告。完成后停止，不生成数据、不启动 fit、不 continuation、不修改包源码。

