# Route-B 含协变量 36-fit 先导：本地与 Linux 启动验收

记录时间：2026-09-27 12:22 +08:00  
协议 ID：`COVARIATE_ROUTEB_PILOT_V1_20260927`

## 本地验收

- runner 与针对性测试脚本解析通过（2/2）；
- runtime manifest 与 `FIT_MANIFEST_CANDIDATE_V1.csv` 的 36 个字符串 `fit_id`、起点编号、拟合种子和单核字段逐项一致；
- 三个臂各 12 条，12 个拟合种子各在 C0/C1/C2 出现一次；
- Route-B 生成器、独立周期 Splinets 均值和缓存后的协变量更新三个相关 testthat 文件全部通过；
- 从零完成本地 `C0 -> C1 -> C2 -> truth-free selection` 微型 smoke：3/3 fit 完成，3/3 臂各冻结一个真值盲 winner；
- 针对性验收 25/25 通过，包括响应分解、C0/C1 配对、C1/C2 原始 Y/time 复用、warning 保留、每 fit 单核和真值隔离；
- 本地微型输出位于 Git 仓库外的临时目录，未纳入正式结果。

## Linux 原生验收

- ECS：`i-bp1egmsairexm9w1rnks`，杭州，8 vCPU、约 16 GiB；
- 上传包 SHA256：`055fa8bfc6ae779a719dc01a3520b6a6f85dc6cbe8e8d78f47d2bdac9eccfaff`，服务器复核一致；
- 独立源码根：`/root/multiFSYNC-covariate-routeb-pilot-20260927-v1`；
- 独立 R 库：`/root/multiFSYNC-covariate-routeb-pilot-lib-20260927`；
- Linux `--action=check` 通过：data=3、fits=36；
- Linux 原生微型 smoke 通过：3/3 fit、3/3 winner、25/25 针对性检查；
- Linux 微型 smoke 峰值常驻内存 116656 KiB，退出状态 0；
- 全尺寸短 sweep 预检：3 sweep 18.210 s，8 sweep 47.305 s，边际 5.819 s/sweep；
- 估算单条 499-sweep 路径约 2903.681 s，36 条、8 路并发约 3.630 小时，不含尾部与 I/O；处于预注册的 3.5--5.5 小时拟合预算内。

## 正式先导启动边界

- 3 个正式观察臂已生成并密封；
- 36 条正式拟合已于本记录形成前以 8 路独立 fit 并发启动，每条内部 `n_cpus=1`；
- 启动时 supervisor PID 为 2629，首批 8 条拟合进程存在，内存使用约 1.6 GiB/14 GiB，未见 `FIT_ERROR.txt`；
- 完成后 runner 只进行臂内最大有效普通 T=1 ELBO 真值盲选择，并冻结 3 个 winner；
- supervisor 结束后调用既有 RAM 身份执行节省停机；
- `truth_unseal_authorized=FALSE`；
- `scientific_evaluation_started=FALSE`；
- 本记录不是科学结果报告，也不对拟合质量作结论。
