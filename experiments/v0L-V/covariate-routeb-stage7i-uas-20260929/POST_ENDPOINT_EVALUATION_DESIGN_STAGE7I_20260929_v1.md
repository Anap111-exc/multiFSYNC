# Stage 7I U/A/S 端点科学评价适配方案（未执行）

日期：2026-09-29。本文只规定后续评价的接口和判定口径；编写时未打开 Stage 7I 真值包、未评价任何 Stage 7I 端点。U、A、S 均为同一 `covrb7c_01__C1__start06` 的 sweep 50 诊断续跑，不属于 Stage 7C 的 12 起点正式候选池，也不据其结果重选 winner。

## 解封前的硬门槛

评价器必须先在**不加载真值**的阶段验证：

1. Stage 7I 的 U 全段 `u_parity_v1/PARITY_PASS.txt` 存在，51–499 状态、400 个普通 ELBO 和端点与冻结 Stage 7G 原路径通过既定门槛。
2. 三条分支各有 499 sweep 完成标记、`fit.rds`、51–499 原始快照、普通 ELBO、盲层轨迹和退出状态。失败或非有限分支也必须在盲层索引中保留；不得只评价成功者。
3. 建立三行 `BLIND_ENDPOINT_MANIFEST.csv`，固定列为 `branch`（顺序 U/A/S）、`source_fit_id`、`data_id`、`fit_path`、`fit_sha256`、`snapshot_path`、`snapshot_sha256`、`blind_metrics_path`、`blind_metrics_sha256`、`status`、`final_ordinary_elbo`、`output_eligible`。路径写成相对 Stage 7I 输出根目录的规范路径；无端点时路径和 ELBO 明确留空，不补跑。绑定观察包、Stage 7G 输入和快照、Stage 7I 源码归档、runner、构造脚本、评价器源和 R 环境的 SHA256。
4. 验证 U/A/S 的数据、维度、预算和停止规则相同，所有盲层表及三行 manifest 的哈希固定在 `BLIND_ENDPOINTS_FROZEN.txt`。marker 必须写明 `truth_used_for_initialization=FALSE`、`truth_used_for_endpoint_selection=FALSE`、`formal_candidate_pool=FALSE`、`diagnostic_case_selected_after_prior_truth_review=TRUE`。即使某一分支异常，盲层仍要先冻结完整状态清单，再决定哪些端点可评价。
5. 核对下表旧输入和评价依赖的 SHA256；任一不符即停。只在上述检查均通过之后读取既有 Stage 7C 执行输出根目录下的 `data/covrb7c_01__C1/sealed_truth/truth_bundle.rds`，并核对其 SHA256 为 `2101ae33f95597b8756437e673a93e690a9d5e7872199307fe6735ed38b46934`、bundle class、data ID、生成器方法和 C1 的 `d=2`。该 hash 来自 Stage 7H 的真值输入绑定，不是在本方案中重新打开真值取得。

Stage 7G 对应的观察包 SHA256 为 `2467e559b055ea9efdac3e8352df1cebce3bcd6940c2260dad20030af552ddd6`，参考 fit SHA256 为 `b2166588d932c4432261e026745f1d1d8185dc0b00f741125d6a779de5d7d761`，快照 SHA256 为 `e85ce3fdd8b2770674232e11cdab5453f61b30271977a41cb08d54b6383621c2`。这些只验证输入身份；Stage 7I 三个新 fit 的哈希须在盲层完成后填入，不得预填。

## 原评价器依赖与适配边界

Stage 7C 的 `evaluate_covariate_routeb_stage7c_20260928_v1.R` SHA256 为 `1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deadbacbbc43b8`。它的前置条件写死为 108 fit、9 个数据 ID、每池 12 个起点、Stage 7C 真值盲 winner 与跨 C0/C1/C2 配对（78–203 行），因此**不能**把三条 Stage 7I 分支装作 Stage 7C manifest 后直接运行。评价科学函数和逐端点判据应保持相同，批量门禁与 winner 汇总需要单独适配。

| 依赖（相对 `experiments/v0L-V`） | SHA256 |
| --- | --- |
| `g12-dimension-overspec-stage6a-20260831/stage6a_common_20260831_v1.R` | `bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6` |
| `g12-reachability-causal-audit-stage6h-20260913/stage6h_evaluation_helpers_20260913_v1.R` | `5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221` |
| `rc1_snapshot/对比实验/v0L阶段A收口与v2冻结准备_20260812/v0lv_candidate2_runtime.R` | `6f0eb7d9c1c01ec3d83074eee43b167cec96c53e7d681ed4ec852f4f791f76bb` |
| 同目录 `r_route_v3_candidate2_runner.R` | `e7d5808bfdd6af756c840eb84b02d87fc3883823742bd93675d1d4d5451003b6` |
| 同目录 `v0lv_candidate2_controller.R` | `ae88153e32bf15455e03e98ade5f64a91a06a47c2947adc38f3da2472a22c5ab` |
| 同目录 `v0lv_candidate2_evaluation.R` | `5ef50d51572e946e158a6b67ee28601c4b8d8847f9830fa11d881e682df24d1e` |
| `rc1_snapshot/对比实验/v0G中性缩减重拟合_20260803/v0g_evaluation.R` | `6b765191657f88db46a027a95acd3917c602b35cb6cee8947d917652f9c6b271` |
| `rc1_snapshot/对比实验/v0H载荷优先结构评估_20260803/v0h_evaluation.R` | `45639c2d9f14e4df0187d519a62837ef3b3be090ed315123cd97c725f8576dd4` |
| `rc1_snapshot/对比实验/v0I真值初始化短诊断_20260805/v0i_truth_diagnostic.R` | `127615551ef1e4dbe14b8c486e3e3ad27ba45eea77f764c7ed7c2fba8318b194` |
| `rc1_snapshot/对比实验/v0J无真值可达性析因_20260808/v0j_evaluation.R` | `3c91af1e07b0ba1aac4539a2c1e28d8abcfe2f5ee239a2a2daed7fcb2714a833` |

评价适配器应加载上述 SHA 绑定的 Stage 6A/6H helper 和 rc1 snapshot，调用**未改动**的 `candidate2_load_frozen_multi_evaluator()` 与 `s6h_evaluate_fit()`。Stage 7C 在调用后者前还执行其 361–465 行的 `install_dense_residual_override(evaluator_environment)`；它修正旧历史 dense residual 的均值装配。遗漏此步骤会使部分 Stage 7C 数值不可比。Stage 7C 的 `mean_recovery()`、`covariate_and_signal_metrics()` 和 `add_recovery_flags()` 位于同一可执行脚本，不能 `source()` 整份脚本以取函数，因为其顶层会立即执行 108-fit 门禁和真值读取。适配器可按已绑定 SHA 解析脚本，仅提取并求值已列名的函数定义（257–788 行），然后显式执行 override；或在 Stage 7I 新文件中逐字复制这些函数并以 U 历史行复核。两种方式均不得修改 Stage 7C 原文件。

`s6h_evaluate_fit()` 固定把 `evaluation_contract$continuation` 写成 `FALSE`（Stage 6H helper 第 180 行），而 Stage 7I 确实从 sweep 50 私有续跑。适配器应保留原函数输出供对照，**另加** `stage7i_private_resume=TRUE` 和 `source_sweep=50`；最终 Stage 7I 契约不得把旧字段单独显示为本次未续跑的证据。`truth_used_for_fit_stopping_or_elbo_selection=FALSE` 不等于诊断案例选择完全真值盲；后者已看过 Stage 7H 的结果，必须单列记录。

## 固定逐端点评价

对 U、A、S 全部可读端点依次运行同一函数和相同 `>=0.5` factor PPI 门槛，保存 `evaluation.rds`、`summary.csv`、`candidate_activity.csv`、`feature_component_errors.csv`、`factor_kernel_errors.csv`、`process_contribution_by_role.csv`、`loading_scale_by_block.csv`、`dense_signal_metrics.csv`、`BETA_COEFFICIENT_ERRORS.csv`、`mean_recovery_by_study.csv` 和 `evaluation_contract.csv`。任何一支评价错误也写出 `EVALUATION_ERROR.txt` 并保留该支；不能只重跑或只汇总成功者。

三个主结构标志使用 Stage 7C 第 772–788 行的**原定义**，不添加事后阈值：

```text
ppi_factor_counts_exact =
    selected_shared == truth_L_f
    AND selected_specific_study_1 == truth_L_s_1
    AND selected_specific_study_2 == truth_L_s_2

core_truth_recovered =
    global_correct == truth_global_direction_count
    AND global_misplaced == 0 AND global_missing == 0
    AND study_correct == truth_study_role_count
    AND study_misplaced == 0 AND study_missing == 0

exact_thresholded_factor_structure =
    core_truth_recovered AND ppi_factor_counts_exact
    AND global_extra == 0 AND global_duplicate == 0
    AND study_candidate_extra == 0 AND study_candidate_duplicate == 0

full_true_feature_component_coverage =
    feature_components_matched == feature_components_total
```

本 C1 数据的真方向数为 3、研究角色数为 4、真特征成分数为 8；这三个计数仍须从所绑定真值和 frozen evaluator 校验，不能按运行结果改写。`core`、`exact`、`8/8 coverage` 分列报告；“完整恢复”至少要求三个标志同时为真，且单列报告 shared、B1、B2 的 `process_contribution_by_role` 与 `loading_scale_by_block`、`R_A/R_B1/R_B2`、跨块 `L_*`、全信号 NRMSE。`core/exact` 本身不保证 shared 过程或完整信号误差达到某个新阈值；本诊断不在看到结果后发明损伤门槛。

`dense_signal_nrmse_ppi_selected` 在历史 helper 中是**不含完整协变量贡献的残差型指标**，不能充作 Stage 7C 的 `obs_full_signal_nrmse_selected` 或 `dense_full_signal_nrmse_selected`。后二者须调用绑定的 Stage 7C `covariate_and_signal_metrics()`，并分别输出。普通 ELBO 只从各支固定的 400 个 T=1 值取最后一项；注入瞬间条件目标不是完整 ELBO。strict convergence 仅报告状态，不参与结构标志。

## 基线复核与交付

先用同一适配器评价 U，核对其 `fit.rds` 与 Stage 7G/Stage 7C 原 winner 的端点复现，然后把 U 的 `summary.csv`、`process_contribution_by_role.csv`、主结构标志和完整信号指标同 Stage 7C 历史的 seed01/C1/start06 行逐字段对照。绑定的 Stage 7C 基线文件 SHA256：`ALL_START_SCIENCE.csv` 为 `fcb9a4db7a8ae47062f8c73435e737b69448695b245e4050e9a8b2ef74554081`，该 fit 的 `summary.csv` 为 `d9aa10a44e6423865b0b035f7c3805df77d460cd3ab1d5122d72f112891da71f`，`process_contribution_by_role.csv` 为 `5daabb2b29316ee3b5d486b653606b9ce5893b098e7ea25c2e47d8a1c6f51bb3`。身份、计数和布尔标志须完全一致；数值用事前固定的 `1e-8 + 1e-8×|参考|` 门槛。若 U 的科学评价与原历史行不一致，应停止 A/S 结构解释并定位包装/后处理差异，不靠调整阈值通过。

最终输出 `STAGE7I_ENDPOINT_SCIENCE.csv`（一行一支，含失败与不可评状态）、`STAGE7I_ROLE_DAMAGE.csv`（一行一支一真实角色）、`STAGE7I_PAIRED_DIFFERENCES.csv`（A−U、S−U、S−A 的普通 ELBO 与连续误差，布尔结构结果保持原值）、`EVALUATION_INPUT_SHA256.csv`、`U_BASELINE_SCIENCE_PARITY.csv` 和完成/错误 marker。真实结构、NRMSE、ISE、R/P/L 只用于**事后评价**，绝不反馈到初值、续跑、选择或提前停止。即使 A 或 S 在 seed01 成功，也只支持本快照的诊断结论，不形成跨 seed 或正式方法成功率。
