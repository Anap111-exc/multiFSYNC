# Stage 6E-B：条件预测 log-score 可行性验证协议 v1

日期：2026-09-10  
状态：development feasibility；不是冻结协议，不授权正式拟合或真值解封。

## 1. 本阶段只回答什么

本阶段只判断：在 `d=0`、整段变量向量按受试者中间时间点留出的既有
Stage 6B 设计下，能否从 multiFSYNC 当前拟合对象保存的变分后验一、二阶矩，
稳定计算一个条件预测 log-score。

以下事项明确不在本阶段改变：

- Stage 6B 已预注册的 held-out RMSE/NRMSE 接口；
- G12、多起点、预得分、退火、普通 ELBO 或 200/400 sweep 规则；
- 模型、先验、CAVI 和 ELBO；
- 任何已冻结或历史实验结果；
- 正式维度选择主准则。

本阶段只运行解析、针对性单元测试和一次微型 smoke，不运行新全尺寸拟合。

## 2. 预测对象与近似边界

在训练时间点拟合后，对同一受试者留出的时间点 `t*` 预测完整的 `p` 维响应
向量。当前均值场变分分解使各参数块独立，但载荷与因子过程的乘积并非高斯，
残差方差也具有后验不确定性。因此，本阶段计算的是：

> 与现有变分后验的一、二阶矩匹配的多元高斯近似 log-score，而不是精确的
> 后验预测密度。

该限定会同时写入机器结果字段
`exact_posterior_predictive_density=FALSE`，避免将近似分数误报为精确 Bayes
posterior predictive score。

## 3. 一、二阶矩公式

以一个共享因子为例。令载荷向量为 `a`，该受试者在 `t*` 的因子过程标量为
`H(t*)`。均值场下 `a` 与 `H` 独立，且不同变量的载荷变分块相互独立。记

\[
\bar a=E_q(a),\qquad v_a=E_q(a^2)-\bar a^2,
\]

\[
\bar H=E_q\{H(t^*)\},\qquad
v_H=\operatorname{Var}_q\{H(t^*)\},\qquad
s_H=E_q\{H(t^*)^2\}.
\]

则该因子贡献的预测均值与协方差是

\[
E_q\{aH\}=\bar a\,\bar H,
\]

\[
\operatorname{Cov}_q(aH)
=v_H\bar a\bar a^\top+s_H\operatorname{diag}(v_a).
\]

因子过程写成 `H=z^T u(t*)`。令

\[
E_q(zz^\top)=\Sigma_z+\mu_z\mu_z^\top,
\]

\[
E_q\{u(t^*)u(t^*)^\top\}
=h(t^*)h(t^*)^\top+
\operatorname{diag}\{v_{u,1}(t^*),\ldots,v_{u,M}(t^*)\},
\]

其中 `h_m(t*)=c(t*)^T mu_nu,m`，
`v_u,m(t*)=c(t*)^T Sigma_nu,m c(t*)`。于是

\[
E_q(H)=\mu_z^\top h,
\qquad
E_q(H^2)=
\sum_{m,m'}E_q(z_mz_{m'})E_q(u_m u_{m'}).
\]

均值函数系数不确定性、每个共享/特异因子的上述贡献以及残差后验均值方差相加，
得到

\[
\Sigma_{\mathrm{pred}}=D+UU^\top,
\]

其中 `D` 为正对角矩阵，`U` 每列对应一个共享或特异因子的
`sqrt(v_H) * mean_loading`。该结构允许通过矩阵行列式引理和 Woodbury 恒等式，
只求解因子数阶的小矩阵，而无需构造或分解 `p x p` 矩阵。

每个留出向量的主实验量候选为

\[
\mathrm{NLPD}(y^*)/p=
\frac{1}{2p}\left\{p\log(2\pi)+\log|\Sigma_{\mathrm{pred}}|
+(y^*-\mu)^\top\Sigma_{\mathrm{pred}}^{-1}(y^*-\mu)\right\}.
\]

越小越好。另输出逐变量边际 NLPD、95% 边际区间覆盖率、标准化平方误差，
它们只用于检查尺度和校准，不在本阶段新增成功门槛。

## 4. 验证门槛

可计算性只有在以下检查全部通过时成立：

1. 新增 R 文件全部可解析；
2. 无低秩项时与独立高斯闭式计算一致；
3. Woodbury/行列式引理结果与显式稠密 Cholesky 结果一致；
4. 双线性因子过程的一、二阶矩与 Monte Carlo 检查一致；
5. 已知高斯分布的 95% 边际覆盖计算回到约 95%；
6. 微型真实 fit 的所有预测方差、joint/marginal NLPD 均有限；
7. 微型 fit 不需要无法解释的方差截断或负方差修补；
8. 新预测均值与原 Stage 6B `factor_mode="all"` 的点预测一致；
9. 整条路径不读取生成真值，每条 fit 仍固定 `n_cpus=1`。

微型 smoke 的经验覆盖率只做描述。由于受试者、变量数和迭代预算都很小，
它不能用于判断模型已经校准，也不能据此把 log-score 升为维度选择主准则。

## 5. 决策规则

- 若上述门槛有一项失败：停止；RMSE/NRMSE + one-SE 保持为当前候选方案，
  并记录 log-score 失败原因。
- 若全部通过：只说明“从现有拟合对象稳定可计算”。下一步仍需在不看真值的
  配对维度候选上进行小规模判别力验证，之后才能比较 log-score 与 RMSE。
- 只有后续验证显示 log-score 在不同拆分上更稳定，并且改善校准/维度排序而
  不牺牲覆盖，才讨论是否修订预注册主准则；不得在本阶段自动替换 RMSE。
- 若后续判别力不佳：保留“重复配对时间块 NRMSE + one-SE”，明确其选择的是
  预测有效维度，不宣称一致恢复生成真值的 `L/M`。

## 6. 与 Jaoua 的维度处理的关系

Jaoua 文中以 `Q` 表示因子数、以 `L` 表示每个因子的 FPCA 截断数；对应
multiFSYNC 的 `L_f/L_s` 与 `M_f/M_s`。她没有为 FPCA 模态使用 MGPS、CUSP
或离散 inclusion indicator。她用较大的 `L_max` 拟合，通过模态特异的样条
方差参数作连续收缩，随后将后处理分数的经验方差作为特征值，按累计 PVE 平台
或高阈值（文中例为 99%）截断。因子数 `Q` 则另用 loading spike-and-slab
推导 factor-PPI 后按阈值保留。

本 Stage 6E-B 的 log-score 是外部候选维度比较工具，不改写上述内部连续收缩，
也不把预测维度选择混同为 posterior rank selection。
