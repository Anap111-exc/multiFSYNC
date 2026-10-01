# Stage 7C ECS 启动记录

日期：2026-09-27  
状态：`FORMAL_STAGE7C_FITS_RUNNING`  
启动时间：2026-09-27T13:02Z（北京时间约21:02）

## 服务器与路径

- ECS：`i-bp1egmsairexm9w1rnks`，`cn-hangzhou`；
- 配置：8 vCPU、约14 GiB可用内存；
- 源码根：`/root/multiFSYNC-stage7c-20260927-v1`；
- 隔离 R 库：`/root/multiFSYNC-stage7c-lib-20260927-v1`；
- 结果根：`/root/multifsync-results/covariate_routeb_stage7c_20260927_v1`；
- supervisor PID：`2494`；
- 外层 worker：8；每条 fit 内部 `n_cpus=1`。

## 启动前验收

- 上传归档：`multiFSYNC-stage7c-20260927-v1.tar.gz`；
- 归档大小：166645 bytes；
- 本地与服务器 SHA256：`93b5fb86fc9cb687154502d00a0920979f51aa7f593966eabb462b2d73fcf924`；
- 13个归档内关键文件逐项 SHA256：13/13 PASS；
- 独立安装 multiFSYNC 0.3.0.9000：PASS；
- Linux registration check：`data=9 fits=108` PASS；
- Linux 微型 smoke：3/3 fit COMPLETE 且 objective-eligible；
- Linux 定向测试：22/22 PASS；
- 微型 truth-free winner：3/3；
- 微型 smoke 未解封真值、未运行科学评价。

## 正式数据与 benchmark

- 3个逻辑数据种子、9个密封观测包生成完成；
- 正式 fit manifest：108行；
- 配对 QC：54/54 PASS，0 FALSE；
- 正式拟合启动前不存在 fit 目录；
- benchmark 的3/8普通 sweep轨迹可复现；
- 估计每个普通 sweep：5.5516秒；
- 仅按 benchmark 线性外推的8-worker下界：约10.39小时；
- 结合36-fit先导的真实退火与长尾运行时间，继续使用15--18小时保守墙钟预算。

## 启动状态

启动核验时：

- supervisor 存活；
- 8个 fit 目录已建立；
- 9个 Stage 7C Rscript 进程（1个批处理监督进程加8个fit子进程）；
- 内存约1.6/14 GiB；
- load average约6.47；
- 无即时错误。

supervisor 结束后会写入 `SUPERVISOR_EXIT_STATUS.txt`，并通过 `multifsync-self-stop` RAM profile 调用 `StopInstance` 的 `StopCharging` 模式。只有108条拟合全部进入 COMPLETE/ERROR 终态时才执行真值盲选择；error记录保留且不选择性重跑。

## 未授权事项

- `truth_unseal_authorized=FALSE`；
- `scientific_evaluation_authorized=FALSE`；
- `continuation_authorized=FALSE`；
- 未启动任何 Stage 7C 真值评价；
- 未修改模型、先验、CAVI、ELBO、初始化或选择规则。
