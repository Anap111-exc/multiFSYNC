# iid-common-Splinets-Stiefel 真值生成器开发协议 v1

日期：2026-09-21  
状态：development candidate；不是冻结模拟协议

## 1. 目的和边界

Stage 7B 表明，旧 `jaoua_splinets` 生成器中的固定真值模板与共享、研究1特异、研究2特异角色绑定，导致连续域误差跟随特定模板移动。本轮只开发一个角色可交换的真特征函数生成选项，不修改模型、先验、CAVI、ELBO、初始化、停止规则、载荷生成规则或评价公式。

本轮只允许解析测试、单元测试、完整生成器微型 smoke 和不拟合的几何审计；不运行任何模型 fit。

## 2. 新方法

公共方法名：`iid_common_splinet_stiefel`。

使用已经确定的 R 包 `Splinets`：

- `Splinets::splinet()` 构造公共连续正交样条系统；
- `Splinets::lincomb()` 构造随机线性组合；
- `Splinets::gramian()` 核验连续内积；
- `Splinets::evspline()` 在任意观测时点直接求值。

仅使用 base R 的 `qr()` 和 `qr.Q()` 取得 Stiefel 系数，不另建样条或正交化框架。

令

\[
\boldsymbol s(t)=(s_1(t),\ldots,s_{J_{\rm truth}}(t))^\top,
\qquad
\int_0^1\boldsymbol s(t)\boldsymbol s(t)^\top dt=I.
\]

v1 固定：

- 三次样条；
- 非周期边界；
- 区间 `[0,1]` 上等距 knots；
- `J_truth=8`；
- `Splinets >= 1.5.1`。

对每个完整函数因子块 `q` 独立生成

\[
Z_q\in\mathbb R^{J_{\rm truth}\times M_q},\qquad
Z_{q,km}\overset{iid}{\sim}N(0,1),
\]

以薄 QR 得到 `Q_q`，并定义

\[
\boldsymbol\phi_q(t)=Q_q^\top\boldsymbol s(t),
\qquad Q_q^\top Q_q=I_{M_q}.
\]

因此每个块内部连续 `L2` 正交。不同块的 `Q_q` 独立生成，不强制块间正交，也不基于拟合质量或真值恢复筛选函数实现。所有共享与特异角色来自同一分布，支持任意总因子块数和不同的 `M_q`，只要求 `M_q <= J_truth`。

## 3. 与载荷生成的关系

函数真值方法与载荷真值方法是两个独立开关。本轮不合并或修改：

- `identified_loadings=TRUE`：强识别载荷场景；
- `sparsity_mode="jaoua_beta_iid_normal"`：Jaoua 式独立 Beta--Bernoulli 随机载荷场景。

严格的 Jaoua 独立支持不保证 Gram 正交。未来若加入“Beta 校准的强识别载荷”，必须另取方法名并明确它不是 Jaoua 原样机制。

## 4. 向后兼容要求

1. `legacy_dense_gram` 继续为默认值。
2. `jaoua_splinets` 的实现和限制保持不变。
3. 新选项不得改变旧方法的 RNG 序列或返回语义。
4. `truth_basis_seed` 必须继续隔离函数生成随机性。
5. 观测点真函数必须由 Splinets 对象直接求值，不得从 dense grid 线性插值。
6. 保存公共 knots、`J_truth`、Stiefel 系数、连续 Gram、网格 Gram 和 Splinets 版本。

## 5. 500次不拟合几何审计

固定 `500` 次独立函数系统抽样，每次生成3个角色块，每块 `M=2`，`J_truth=8`，评价网格201点。每次只计算指标并丢弃函数对象，不创建500个代码或RDS文件。

审计指标：

- 块内连续 Gram 最大误差；
- 201点梯形积分 Gram误差；
- 块间最大典型相关；
- 二阶差分粗糙度；
- 两端各10%区间的边界能量比例；
- 投影到当前 `K=5` O'Sullivan 拟合空间的相对 `L2` 误差下界；
- 三个角色块的分布可交换性。

预先规定的有效性门槛：

1. 所有连续 Gram 最大误差不超过 `1e-8`；
2. 全部指标有限；
3. 不出现典型相关达到 `1-1e-8` 的数值重复子空间；
4. 块间典型相关中位数大于 `0.1`，确认实现没有误加全局正交；
5. 各角色块均值范围除以合并标准差不超过 `0.25`；
6. 投影误差下界总体中位数不超过 `0.25`，95%分位数不超过 `0.60`。

第5项分别检查粗糙度、边界能量和投影误差。个别高重叠或困难实现不作选择性删除；只报告尾部比例。

## 6. 后续决策

若单元测试和几何审计通过，才进入独立的小规模拟合桥接；若失败，只调整全局生成超参数并新建协议版本，不按模型恢复结果挑选函数种子。本协议本身不授权任何拟合。
