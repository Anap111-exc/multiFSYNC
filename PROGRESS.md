# multiFSYNC — 多研究函数型因子模型实现进度

## 项目信息
- **目标**：多研究函数型因子模型 R 包，支持共享因子、研究特异因子、协变量效应
- **基于**：bayesSYNC (https://github.com/hruffieux/bayesSYNC)，GPL-3
- **工作目录**：`D:\文档\Factor model\multiFSYNC\Rcode\multiFSYNC\`
- **索引约定**：研究维度永远最外层 `[[s]]`，个体 `[[i]]`，因子 `[[l]]`，FPCA 分量 `[[m]]`

> **2026-07-24 当前说明**：下列早期阶段记录保留为开发历史，其中旧 ELBO 公式与“实验性诊断”结论已被 Phase 11 取代。当前实现只在完整 `T=1` sweep 后计算普通 ELBO；仅当 `lambda_orth=0` 且该阶段未触发自适应 jitter 时允许 ELBO 停止。`T>1` 的 PPI 路径仍采用论文附录公式，并以最终 `T=1` 结果为准。

## 当前状态

| 阶段 | 状态 | 完成时间 |
|------|------|----------|
| Phase 0: 环境准备 | 完成 | 2026-05-15 |
| Phase 1: R 包骨架 | 完成 | 2026-05-15 |
| Phase 1.1: generate_data.R | 完成 | 2026-05-15 |
| Phase 2: 初始化与预计算 | 完成 | 2026-05-15 |
| Phase 3: 模块实现 | 完成 | 2026-05-15 |
| Phase 4: 主循环集成 | 完成 | 2026-05-15 |
| Phase 5: 后处理 | 完成 | 2026-05-15 |
| Phase 6: 测试套件 | 完成 | 2026-05-15 |
| Phase 7: 端到端验证 | 完成 | 2026-05-15 |
| Phase 10: 旧 ELBO 诊断实现 | 已由 Phase 11 取代 | 2026-07-23复核 |
| Phase 11: 完整普通 ELBO 修复 | 完成并通过小样本/回归测试 | 2026-07-24 |

## Phase 11：完整普通 ELBO 修复

- 数据项使用最后一次载荷更新后的直接完整后验期望 RSS，并由误差方差更新缓存供 ELBO 复用。
- O'Sullivan 线性/非线性两块、标准正态计算得分、共享/特异 spike-and-slab、两种 omega 层级及全部 Half-Cauchy IG--IG 层级均纳入命名分量。
- `update_beta`、`update_eigenfunctions`、`update_zeta`、`update_xi` 与 `update_loadings` 的耦合坐标改为工作副本上的 Gauss--Seidel 顺序更新。
- 高斯精度首先对未扰动矩阵作 Cholesky；失败时才逐级增加 jitter，并输出位置与次数。`T=1` 触发 jitter 或 `lambda_orth>0` 时自动改用参数停止。
- 新增 `test_elbo_components.R` 与 `test_elbo_monotonicity.R`。标准
  `testthat` 全套回归通过；高耦合场景覆盖
  `d=2, L_f=L_s=2, M=2`，8 个 ELBO 点无下降，并逐元素核对最终
  RSS 缓存、两种 omega 的 ELBO 装配以及 T=1 jitter 自动降级。

## 已创建文件

### 包结构
- `DESCRIPTION`、`NAMESPACE`
- `R/utils_multi.R` — 工具函数集 (~890行)

### 核心
- `R/generate_data.R` — simulate_multi_study_data()
- `R/multi_core.R` — bayesSYNC_multi() + core + compute_elbo_multi()
- `R/orthonormalise_multi.R` — SVD正交化后处理

### CAVI 更新模块 (8个)
- `R/update_variance.R` — 5方差 + 5辅助变量
- `R/update_omega.R` — 共享/特异 omega
- `R/update_mu.R` — 均值函数 nu_mu
- `R/update_beta.R` — 回归系数 nu_beta (跨研究聚合)
- `R/update_eigenfunctions.R` — nu_phi + nu_psi (跨因子交叉项)
- `R/update_zeta.R` — 共享得分 zeta (H矩阵)
- `R/update_xi.R` — 特异得分 xi
- `R/update_loadings.R` — a+b载荷+PPI (熵项不乘c)

### 测试文件 (9个)
- `tests/testthat/test_degenerate_jaoua.R`
- `tests/testthat/test_degenerate_T.R`
- `tests/testthat/test_dimensions.R`
- `tests/testthat/test_ppi_temperature.R`
- `tests/testthat/test_gaps.R`
- `tests/testthat/test_elbo_components.R`
- `tests/testthat/test_elbo_monotonicity.R`

---

## Phase 6 测试套件：逻辑与结果

### test_degenerate_jaoua.R — JAOUA 退化测试

**测试目标**：验证 S=1, K_s=0, d=0 时多研究版本能正常运行。

**测试逻辑**：
1. 同种子生成单研究数据 (N=10, p=3, Q=1, L=2)
2. 分别用 bayesSYNC 和 bayesSYNC_multi 拟合 (maxit=20, anneal=NULL)
3. 验证两者 ELBO 有限、PPI 在 [0,1]

**结果**：通过。bayesSYNC ELBO=-498, multiFSYNC ELBO=-54。ELBO 不精确匹配已知原因：PPI公式在c=1时等价，但 sigma^2 IG 熵项和辅助变量处理方式不同导致 ELBO 量级差异。两者 PPI 均在合法范围。

### test_degenerate_T.R — 温度/T 退化测试

**测试逻辑**（3个子测试）：

| 子测试 | 参数 | 标准 | 结果 |
|--------|------|------|------|
| T=1 标准CAVI | anneal=NULL | ELBO有限, PPI in [0,1] | 通过 |
| T→0 退火 | anneal=c(1,1.9,20) | 退火收敛, PPI合法 | 通过 |
| T>=2 越界 | anneal=c(1,2.0,10) | 抛出错误(kappa_a=2c-1<=0) | 通过 |

### test_dimensions.R — 维度检查

**测试逻辑**：运行 S=2, K_f=2, K_s=1, d=1, L_f=2 完整模型 (maxit=5)，逐项断言维度。

| 参数 | 预期维度 | 结果 |
|------|---------|------|
| mu_q_nu_mu | list[[s]][[j]], (K+2)x1=8 | 通过 |
| mu_q_nu_beta | list[[j]][[r]], (K+2)x1 | 通过 |
| mu_q_nu_phi | list[[l]][,m], (K+2)xL_f=8x2 | 通过 |
| mu_q_nu_psi | list[[s]][[l]][,m], (K+2)xL_s=8x1 | 通过 |
| mu_q_zeta | list[[s]][[l]], n_s[s]xL_f=5x2 | 通过 |
| mu_q_xi | list[[s]][[l]], n_s[s]xL_s=5x1 | 通过 |
| loadings a | p x K_f = 3x2 | 通过 |
| list_cp_C[[s]][[i]] | (K+2)x(K+2) = 8x8 | 通过 |
| sum_list_cp_C[[s]] | (K+2)x(K+2) = 8x8 | 通过 |

### test_ppi_temperature.R — PPI 温度效应（核心防守）

**测试目标**：严格验证 PPI logit 公式中高斯熵项不被温度 c 缩放。

**子测试 1 — 孤立公式验证**：
```
设定：mu=2, sigma^2=1, E[log omega]=0, E[log(1-omega)]=0
熵项 = 0.5*(4/1 + log 1) = 2

c=1:   logit = 1*(0-0) + 2 = 2  -> PPI = 0.8808
c=0.5: logit = 0.5*(0-0) + 2 = 2 -> PPI = 0.8808 (不变!)
       若错误乘c: logit = 0.5*(0+2) = 1 -> PPI = 0.731
```
结果：通过。该测试验证本文采用的参考论文附录公式；本地 bayesSYNC 代码在高温阶段使用不同实现，不能仅凭此测试将其称为错误。两种写法在最终 `T=1` 阶段一致，正式实验仍需报告敏感性分析。

**子测试 2 — 集成测试**：S=1, K_f=1, maxit=10, PPI in [0,1] → 通过

**子测试 3 — 多因子PPI范围**：S=1, K_f=2, maxit=10, 所有PPI in [0,1] → 通过

### test_gaps.R — 缺口测试

**子测试 1 — orthonormalise 数值验证**：正交化后 integral(phi^2) = 1 (tolerance=0.1), PVE sum=100% → 通过

**子测试 2 — factor_ppi 公式验证**：PPI=(0.2,0.5,0.8) -> factor_ppi = 1 - 0.8*0.5*0.2 = 0.92 → 通过

**子测试 3 — JAOUA T>1 退火**：anneal=c(1,1.9,30), 不与bayesSYNC比绝对值, 自身收敛即可 → 通过

**子测试 4 — 方差迹随c单调递减**：c in {1,5,20}, tr(Sigma)单调递减 → 通过

---

## Phase 4 关键修复：ELBO 从递减到递增

> **已作废的历史诊断**：本节所称 ELBO 仍缺少 Phase 11 后补齐的
> 数据、O'Sullivan、IG--IG、spike-and-slab 等项；这里的升降现象只
> 记录旧调试过程，不能作为当前普通 ELBO 正确性或单调性的证据。

**问题**：ELBO 从 -870 单调递减至 -911。

**根因**（3项）：
1. E[log sigma^2] 计算错误：代码用 `log(E[1/sigma^2])` 近似 `E[log sigma^2]`，但逆Gamma后验期望是 `log(lambda) - digamma(kappa)`，两者差异巨大
2. 缺少 mu 先验+熵正贡献项 (elbo_mu)
3. 主循环未提取 mu_q_log_sigsq_* 传递给ELBO

**修复**：提取 mu_q_log_sigsq_eps/mu, lambda_q_sigsq_eps/mu 传入 compute_elbo_multi, 添加 elbo_mu + elbo_zeta 项

**结果**：ELBO 从 -41.7 上升至 208.6 (迭代3后27/27步全部上升)

---

## Phase 10：ELBO 公式补全与模型对比

> **已作废的历史记录**：本节使用的是 Phase 11 审计前的不完整
> ELBO，只保留用于说明开发过程。以下旧公式、数值差异归因和模型间
> ELBO 结论均不得作为当前实现或论文结论引用；当前依据见 Phase 11
> 与 `理论推导/derivations.md` 第七节。

### ELBO 公式修复

**缺失项**：
- sigma^2 逆Gamma后验熵: (kappa-0.5)*E[log sigma^2] - kappa*log(lambda) + lgamma(kappa) - lgamma(0.5)
- 辅助变量抵消项: bayesSYNC在elbo_y用 -E[1/sigma^2]*(lambda-E[1/a])，同时elbo_mu补回相同量，净效果0

**修复** (compute_elbo_multi):
- elbo_y: 辅助项改为 (lambda - E[1/a]) 对齐 bayesSYNC
- 新增 elbo_sigma 段: sigma^2_eps + sigma^2_mu 的IG熵 + 辅助变量抵消
- 总ELBO = elbo_y + elbo_mu + elbo_loadings + elbo_zeta + elbo_sigma

### multiFSYNC vs bayesSYNC 对比结果

**T=1 (c=1) PPI 对比** (S=1, N=25, p=8, Q=2, L=2, 有结构的数据):
- PPI 相关矩阵: bayes factor 2 ↔ multi factor 1 corr=0.997, 次优 corr=0.503
- 平均 |corr| = 0.75
- ELBO delta = 4297 (bayesSYNC=-74 vs multiFSYNC=-4371)

**极小数据诊断** (N=5, p=2, 3次迭代):
- ELBO delta = 121 (bayesSYNC=-245 vs multiFSYNC=-124)
- 差距在该旧实验中缩小；这不能证明差异来自数据规模或
  `kappa_q` 量级。

**历史观察（相关结论已撤销）**：
- T=1时PPI公式等价，因子结构有部分一致性(corr=0.997匹配)
- 旧 ELBO 的绝对值差异来源当时未被逐项识别。
- 旧实验中的样本量趋势不构成对差异来源的证明。
- “两模型 ELBO 公式本质正确”的旧判断已由 Phase 11 审计撤销。

---

## 关键设计决策

1. **PPI 公式**：高斯熵项 `0.5*(mu^2/sigma^2 + log sigma^2)` 不乘 c，仅先验项 `c*(E[log omega]-E[log(1-omega)])` 受温度调控。c=1 时与 bayesSYNC 等价。

2. **索引约定**：研究维度永远最外层 `[[s]]`。共享参数 (nu_phi, nu_beta) 无研究索引。

3. **先验精度**：blkdiag(inv_Sigma_beta, mu_q_recip_sigsq * diag(K)) 用于所有样条系数。

4. **临界约束**：kappa_q_a = 2c - 1 > 0 => T < 2。

## 已知限制

1. **历史 ELBO 对比不再有效**：Phase 10 的对比基于不完整
   ELBO，不能用于验证当前公式。当前完整普通 ELBO 由分项解析测试、
   当前 RSS 一致性测试和 $T=1$ 单调性测试验证；与 bayesSYNC 的数值
   ELBO 不作跨模型等值主张。

2. **计算效率**：纯R for循环实现，未做 C++/Rcpp 加速。大规模数据建议并行化。

3. **PPI 对比仅在 T=1 完成**：T>1 时两模型 PPI 公式不同，不应做精确对比（各自收敛即可）。

## 13 组公式映射

全部 13 组公式在 formula_code_map.md 中均有映射，Phase 1 验证完成。
