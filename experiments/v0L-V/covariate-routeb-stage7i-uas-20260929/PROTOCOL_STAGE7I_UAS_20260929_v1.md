# Stage 7I：同一冻结快照的 U／A／S 幅值对照与谱重启

日期：2026-09-29。此协议在任何 Stage 7I 续跑前固定。Stage 7H 的[原报告](../covariate-routeb-stage7h-residual-driver-20260928/REPORT_STAGE7H_READONLY_RESIDUAL_DRIVER_20260929_v1.md)和[勘误](../covariate-routeb-stage7h-residual-driver-20260928/ERRATUM_STAGE7H_INTERPRETATION_20260929_v1.md)提供只读归因。Stage 7I 只检验固定 seed01/C1/start06 的同一个 sweep 50 快照，不产生正式候选，不修改 Stage 7C–7H 冻结输入或项目当前工作树的模型文件。

## 问题、输入与边界

Stage 7H 表明，sweep 50 的研究 2 留出残差首变量方向对真 B2 的投影平方为 `0.993548`，但原 B2 的得分—时间函数均值乘积和载荷更新驱动极小。三个 winner 的 `‖h‖₂/sqrt(U_mean)` 约为 `83.893/79.945/86.495`；这是幅值对照的动机，不是纯幅值致因的证明。本实验只用 seed01 的已冻结真值盲 winner `covrb7c_01__C1__start06`，在**同一状态**比较未干预 U、保持原方向的幅值恢复 A、留出残差谱重启 S。

必需输入：Stage 7G 的 `FIT_MANIFEST.csv`、`analysis_v1/INPUT_BINDING.csv`、seed01 `fit.rds` 与 `RAW_STATE_SNAPSHOTS.rds`；Stage 7C seed01 `observation_bundle.rds`；Stage 7G 冻结源包及 runner；Stage 7H 的真值盲驱动表和哈希清单。按已有 manifest 验证每个 SHA256。真值 bundle 在构造、续跑、盲层选择和冻结时**不打开**。所有拟合变量使用冻结的原变量顺序、原观测时间、原基和工作尺度。若任何输入绑定不符，停止。

## 三个固定分支及真值盲初值

记 sweep 50 末研究 2 的 B2 有效载荷为 `b_U=ρ⊙m_U`，时间系数矩阵为 `Ψ_U`，个体得分均值为 `ξ_Ui`，`g_Ui=C_i Ψ_U ξ_Ui`。`ρ` 是保存的 PPI，`m` 是 slab 条件均值。所有分支沿用相同 PPI、slab 方差、函数/得分协方差、方差与 omega 状态、其他块状态、先验、退火表和 CAVI 顺序。

**U：** 不改变 sweep 50 快照任何原始状态。

**S：** 对每名研究 2 受试者构造 `R_i=Y_i−C_i μ_i−Σ_r Z_ir C_i β_r−(C_i Φ ζ_i)aᵀ`；本研究的 B1 属于研究 1，不在此残差中。令 `G_R=Σ_i R_iᵀR_i`，取最大特征值 `λ₁` 及单位首特征向量 `v`，将其最大绝对坐标定为正号。固定目标有效载荷范数 `B*=sqrt(λ₁/N_obs)`，其中 `N_obs=Σ_i n_i`；这是工作尺度规范，不声称它是唯一真实幅值。令 `q_i=R_i v`，`σ_v²=Σ_j v_j²/τ_j`，其中 `τ_j` 是 sweep 50 保存的研究 2 噪声精度。

使用原 O'Sullivan `K+2` 维基、原观测设计与**固定模型先验形状** `P=blockdiag(Σ_beta^{-1}, I_K)`，其中 `Σ_beta` 来自冻结默认超参数；平滑系数固定为 `λ=σ_v²/B*²`。先解跨受试者共同时间均值 `δ=(Σ_i C_iᵀC_i+λP)^{-1}Σ_i C_iᵀq_i`，取 `r_i=q_i−C_iδ`。构造 `D=Σ_i C_iᵀr_i r_iᵀC_i`、`H=Σ_i C_iᵀC_i+λP`，取 `D u=d H u` 的前 `M_s,2=2` 个广义特征方向，按原 51 点网格范数规范并对各列最大绝对网格值定正号，作为初始 `Ψ_S`。各人的初始二维得分为

`ξ_Si=(B*² Ψ_Sᵀ C_iᵀC_i Ψ_S+σ_v² I₂)^{-1} B* Ψ_Sᵀ C_iᵀ r_i`。

随后固定执行**恰好三轮**交替更新：依次按 `m=1,2` 解

`(Σ_i ξ_im² C_iᵀC_i+λP) ψ_m=(1/B*)Σ_i ξ_im C_iᵀ[r_i−B*Σ_{k≠m}C_iψ_k ξ_ik]`，

再用上式更新全部 `ξ_i`。最后逐列中心化得分，按其经验 RMS 规范为 1 并将逆尺度转到对应函数列；同时翻转每列函数与得分的符号，使最大绝对网格函数值为正。最后设 `b_S=B* × v`。不得根据真值、续跑结果或端点反复选择谱秩、惩罚、轮数、符号或尺度。

**A：** 保留原有效载荷方向、每人的原 `g_Ui` 方向及有符号累加模式。设 `b_A=B*b_U/‖b_U‖₂`。先将原 `ξ_U` 的两列各规范为经验 RMS 1、把相反尺度转到对应 `Ψ_U` 列；然后对全部函数列乘同一个正数 `G_S/G_U`，其中 `G_X²=Σ_i‖C_i Ψ_X ξ_Xi‖²`。因此 `b_A` 与 `b_S` 的载荷范数同为 `B*`，`G_A=G_S`，且初始完整后验**均值**产品在观测点的 Frobenius 范数相同。A 的每个 `g_Ai` 与原 `g_Ui` 成同一正比例；A 不改变原方向与跨人有符号相消结构。

在 A/S 中仅替换 B2 的载荷、函数和得分**均值**。保留 `ρ_j` 与 slab 方差 `V_j`，由新有效载荷 `b_j` 设置 `m_j=b_j/ρ_j`，重算 `E[b_j²]=ρ_j(V_j+m_j²)`，并验证 `ρ⊙m=b`；其他二阶矩和协方差保持原值。匹配 `‖b‖` 与 `G` 不保证噪声加权二阶矩、得分精度或 PPI logit 相同，必须逐项记录。`ρ` 不严格正、原 `b_U/G_U` 为零、谱秩不足、广义特征值不稳定、非有限状态或幅值匹配不通过时，停止，不能临时改阈值/惩罚。

## 精确退火续跑门槛

现有 `continuation_state` 明确只支持 `anneal=NULL` 的 T=1；不能用它从 sweep 50 续跑。只在隔离的 Stage 7I 源码复制中加入私有恢复钩子，不改变公开接口、CAVI 更新、先验或 ELBO。原 `anneal=c(1,1.9,100)`、`maxit=499`、`pre_score_sweeps=1`、`n_cpus=1`、固定 C1 维度和 practical control 均保持。恢复 sweep 50 末原始状态及 omega，使用原退火表的 `ladder[51]` 从**全局 sweep 51**运行到 499：剩余 49 次退火更新和 400 次普通更新。sweep 50 尚无普通 ELBO/practical 历史，相关缓冲从空开始；不能把局部 449 步重编号为 1–449。

先以 U 做两步微型 smoke，比对原快照 51、52 的温度及每项原始状态。通过后，U 全程续跑，并与 Stage 7G 的 51–499 快照逐项核对：有限数值要求 `|差|≤1e−8+1e−10×|参考|`，PPI 零/非零及形状精确一致；普通 ELBO 路径逐点绝对差 `≤1e−6`，最终结构/状态与原 winner 相符。原始全程求解器调用计数未保存在快照中，不作端点逐字节门槛；续跑段的 jitter、数值修正仍须记录。**任何关键复现失败，停止 A/S，定位首个分歧，不放宽阈值掩盖。**

快照未保存完整逆 Gamma 的形状、率、对数期望，下一次方差更新会重建它们。因此注入瞬间只比较同一固定方差状态下的**条件退火目标差**，明确标记缺少方差-only 项；sweep 51 更新后再记录完整退火目标。不得把注入瞬间的条件值称作完整 ELBO。普通 ELBO 从 sweep 100 开始。

## 记录、判断与停止

三分支均保存：输入/源码/脚本哈希，快照 50 的完整基线与注入状态，51–499 每个 sweep 的原始状态，51–70 每个块后的载荷/PPI/函数/得分摘要和更新驱动分项，全部温度、条件/完整退火目标、普通 ELBO、完整后验乘积二阶 RMS、其他块均值贡献、PPI logit 分项、得分精度矩阵和 `Q_i^{-1}h_i`、个体间有符号驱动与相消、数值修正及停止原因。若单个诊断的完整分块 trace 无法在预算内保存，至少保存全部 sweep 原始状态，保证之后可只读重建。

固定比较 U/A/S 的同预算端点与早期轨迹。**先锁定真值盲表和哈希，再单独解封真值评价**；真值绝不用于初值、分支选择或续跑端点选择。结构成功须同时评价 B2 支持、载荷方向、完整重构、shared/其他块损伤及跨块分配，不能只看 PPI、ELBO 或块能量。A、S 都成功表示谱方向对本路径恢复不必要；仅 S 成功支持残差引导的组织在该幅值下有额外作用；均短暂增长后关闭只说明两种重置不足，仍须区分初值质量、精度、PPI 与块竞争。高 ELBO 但结构错只表示目标改善，非结构可达性修复。

硬上限是 U/A/S 各 `449` 个续跑 sweep，共 `1347`；不追加起点、谱惩罚网格、迭代或第二数据。输入哈希错误、U 复现失败、状态不合法、非有限值、资源上限或完成固定预算即停止。这个诊断不进入正式候选池，不证明从迭代 0 的真值盲可达性、总体失败率降低、一般识别、跨 L/M 选择或 genuine `d_fit=0` 性能。原冻结文件均保留；放弃本诊断只需忽略 Stage 7I 新目录。

## Fixed implementation gates before any real-data continuation

The isolated initializer uses `P = blockdiag(Sigma_beta^{-1}, I_K)` only as a fixed spectral initialization penalty shape. It is not the current full variational prior precision for `q(psi)`. The 51-point grid norm uses trapezoidal weights on the frozen `time_g`; its Gram is `C_g' diag(w) C_g`. The two selected generalized eigenvectors are normalized on this grid, with deterministic largest-absolute-grid-value sign anchors. No subject is interpolated with an 11-column basis.

Numerical gates are fixed in `build_stage7i_uas_states_20260929_v1.R` before looking at real-data branch outcomes: `rcond(H)>1e-12`, each ALS linear system has `rcond>1e-12`, `d2/d1>1e-8`, `(d2-d3)/d2>1e-8`, initial and final two-column grid Gram `rcond>1e-10`, each required old/new mean and score RMS exceeds `sqrt(.Machine$double.eps)`, and A/S amplitude matching uses absolute and relative tolerances of `1e-10`. A failed gate stops Stage 7I; these constants are not adjusted after inspection. Save the exact gate values and blind diagnostics, including removed pooled-mean energy and the final projected residual SSE.

A and S match effective loading norm, observed-time score-function product norm, and hence initial mean-product RMS. They do not match posterior second moments, precision in the actual driver direction, PPI logit, time organization, or residual alignment. If only S recovers, the result supports an added role for the bundle of residual-guided direction and time organization at this amplitude; it does not isolate which component of that bundle caused recovery. If A also recovers, residual guidance is not necessary for this particular continuation, but this does not identify why the original U path lost mean amplitude.

## Locked scientific readout after all endpoints are frozen

This is a diagnostic on the seed01/start06 sweep-50 state, chosen after earlier truth-assisted analyses. Only the new initializer and continuation branches are truth-blind. Do not describe the whole study as a prospective blind method validation.

After U/A/S endpoints, source hashes, blind tables, and branch-completion markers are frozen, use the Stage 7C evaluator logic in `evaluate_covariate_routeb_stage7c_20260928_v1.R` (SHA256 `1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deacbacbbc43b8`), including its bound Stage 6A/6H support hashes (`bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6` and `5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221`). Preserve the original >=0.5 PPI threshold, matching rules, and no after-the-fact threshold changes. Report `core_truth_recovered`, `exact_thresholded_factor_structure`, `full_true_feature_component_coverage`, matched shared and study-specific roles, all missing/misplaced/extra/duplicate counts, 8/8 component coverage, observed and dense full reconstruction errors, beta errors, and complete contribution errors. Call a branch a complete structural recovery only if exact structure and full component coverage both hold; report core separately. Compare shared and study-1 roles and complete contribution errors against U even when B2 improves. Ordinary ELBO, PPI, or B2 support alone are insufficient.

Record `sum E[b_B2^2]`, tau-weighted B2 loading norm, mean-function and score penalties, precision applied to the actual score driver, signed per-subject driver contributions, and PPI logit components. These differ between A and S despite the three matched raw mean amplitudes. The residual leading eigenvalue contains noise and other residual structure, so `B_target` is a fixed observed-residual amplitude convention, not an estimate of a true loading norm. If only S succeeds, infer only that this fixed S bundle outperformed A at the matched amplitude on this state; do not assign the effect to spectral direction, time organization, or scale allocation individually.

## Injection-time algebraic checks

Let `kappa = G_S/G_U` and keep the same leave-B2-out residual and noise precision at sweep 50. Before any CAVI update, A must satisfy `g_Ai = kappa g_Ui`, `h_A = kappa h_U`, `U_mean,A = kappa^2 U_mean,U`, and `||h_A||/sqrt(U_mean,A) = ||h_U||/sqrt(U_mean,U)`. The initializer asserts these identities to its fixed numerical tolerance and saves the signed subject contributions `c_i = sum_j v_j tau_j <g_i, R_ij>`, their sum, absolute-sum cancellation ratio, total and v-directed driver norms for U/A/S. This is a fixed-state algebraic contrast, not proof about the subsequent optimization path. The Stage 7H seed01 sweep-50 number 83.8931320186653 used sweep-49 noise precision with sweep-50 function/score state. The injection calculation uses the complete sweep-50 state, including its noise precision, and need not equal the earlier coordinate-update number; it must satisfy the U/A identity on the same fixed sweep-50 residual and precision.
