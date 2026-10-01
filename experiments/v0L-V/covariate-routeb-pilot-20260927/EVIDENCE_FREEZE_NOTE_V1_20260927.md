# Route-B 含协变量 36-fit 先导证据冻结说明 v1

证据 ID：`COVARIATE_ROUTEB_PILOT_EVIDENCE_V1_20260927`  
冻结日期：2026-09-27  
状态：`frozen_development_evidence`  

## 冻结对象

本记录固化一个逻辑数据种子 `covrbp_01` 下的三臂、36 条固定 horizon 拟合及其科学评价：

- C0：生成 `d=2`、真实协变量效应为零，拟合包含协变量；
- C1：生成 `d=2`、真实跨研究共享 B 样条协变量效应非零，拟合包含协变量；
- C2：精确复用 C1 的响应和观测时点，但拟合不提供协变量；
- 每臂 12 个配对拟合种子；每条拟合为 random + G12 + 一次预得分 + Jaoua 退火 + 400 个普通 `T=1` sweep；
- 所有 36 条拟合完整且 objective-eligible；无选择性重跑；
- 三个 truth-free winner 在解封真值前按臂内最大有效普通 `T=1` ELBO 冻结；
- 解封后完成 36/36 科学评价，未启动新拟合或 continuation。

精确文件、字节数和 SHA256 见 `EVIDENCE_MANIFEST_V1.csv`；机器冻结状态见 `FROZEN_EVIDENCE_STATUS_V1_20260927.txt`。

## 主要冻结结果

| 项目 | C0 | C1 | C2 |
|---|---:|---:|---:|
| core / exact / full-feature 命中 | 6/12 | 6/12 | 0/12 |
| strict practical convergence | 7/12 | 6/12 | 8/12 |
| ELBO winner | start06 | start06 | start08 |
| winner 是否 core/exact/full-feature | 是 | 是 | 否 |
| winner observed full-signal NRMSE | 0.1450 | 0.1458 | 0.6895 |
| winner dense full-signal NRMSE | 0.2292 | 0.2233 | 0.7247 |

C1 winner 的 beta 总能量比为 0.9784、relative ISE 为 0.0533、平均连续 `L2` cosine 为 0.9764。C0 winner 的伪 beta 能量为 C1 真实 beta 总能量的 1.97%。C2 的 12 个起点虽均保留 `1/1/1` 个因子，但均未恢复正确共享—特异身份，说明遗漏协变量可由错误因子方向吸收，而不是简单表现为因子数异常。

严格收敛与科学正确性继续分离：C0/C1 均同时存在严格但非 core、以及未严格但 core 的起点；C2 有 8 条严格收敛轨迹却无一 core。严格收敛不能替代载荷身份和连续科学质量评价。

## 证据边界

本证据支持：

- 当前 Route-B 生成机制和含协变量更新能够在该逻辑数据种子上恢复非零共享 beta；
- 零效应臂未出现实质性伪 beta；
- 在候选池含正确盆地时，普通 ELBO 在 C0/C1 选中了科学质量优良的端点；
- 遗漏协变量会严重污染残差共享—特异结构。

本证据不支持：

- 将 12 个起点当作 12 个独立数据重复；
- 对不同数据种子、不同协变量分布或不同 beta 形状作普遍推广；
- 声称严格 practical convergence 是科学正确性的必要或充分条件；
- 声称已经完成 `L_f/L_s` 或 `M_f/M_s` 的维度选择；
- 声称本先导是论文最终 Monte Carlo。

本地固化包含协议、源码绑定、runner/test、冻结 winner 哈希、完整科学评价归档、汇总表和报告。36 个完整拟合对象未在本地重复复制；其服务器原始目录为 `/root/multifsync-results/covariate_routeb_pilot_20260927_v1`。因此本地证据足以审计当前报告与评价，但若未来删除服务器实例，仍应先另行备份完整 fit RDS；本记录本身不冒充完整原始拟合归档。

## 后续关系

Stage 7C 将使用三个全新逻辑数据种子独立确认上述结论。Stage 7C 不追溯修改本证据，不在本轮启动，也不得根据本先导的真值结果调整其候选 seeds、起点或选择规则。
