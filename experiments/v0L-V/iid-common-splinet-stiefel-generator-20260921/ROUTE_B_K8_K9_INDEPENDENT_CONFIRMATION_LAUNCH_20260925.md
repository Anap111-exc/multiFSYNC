# Route-B K=8/9 独立确认启动记录（2026-09-25）

- 服务器：Alibaba Cloud ECS `i-bp1egmsairexm9w1rnks`，8 vCPU / 16 GiB。
- 服务器源码根：`/root/multiFSYNC-routeb-bridge-20260922-v1`。
- 输出根：`/root/multifsync-results/route_b_k8_k9_independent_confirmation_20260925_v1`。
- runner SHA256：`b4310d1eb1ae01ad748cadf0e9e7001fa244e445f9cf3a550fa4b5bff9b688a4`。
- 目标包：`multiFSYNC 0.3.0.9000`；R 库显式绑定为 `/root/rb22lib:/root/multiFSYNC-jaoua-paired-lib-20260909`。
- Linux 解析：PASS。
- 运行级检查：PASS（data=3，fits=48，K=8/9，starts=8）。
- 微型 smoke：2/2 COMPLETE、2/2 eligible、2 个真值盲 winner 冻结；最大常驻内存约 118 MiB。
- 全尺寸数据生成：3/3 完成并密封；`DATA_MANIFEST.csv` 3 行、`FIT_MANIFEST.csv` 48 行；未发现 `truth_unseal_authorized=TRUE`。
- K=9 benchmark：3/8 sweep 轨迹前缀一致；估计边际 2.7318 秒/T1 sweep；8 并发投影时间 2.271947 小时。
- 主批量：2026-09-25 约 04:56:44 UTC 启动，supervisor PID 2674；每条 fit `n_cpus=1`，外层 8 并发。
- 启动后检查：8 条 fit 同时运行，0 COMPLETE、0 ERROR（检查时点 2026-09-25T04:59:32Z）。
- 自动停机：RAM role `multifsync-ecs-self-stop` dry-run 通过；批量终止并写盘后等待 60 秒，以 `StopCharging` 停止实例。

本轮到 truth-free selection 冻结即停止。真值仍密封，科学评价、continuation 和跨 K 选择均未启动。
