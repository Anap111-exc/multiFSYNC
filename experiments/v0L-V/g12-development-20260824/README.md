# v0L-V G12 协议提升与 convergence 语义审计

本目录记录 2026-08-24 用户批准的两项 post-v0L-V 开发工作：

1. 将 `gram_unit_energy` 提升为下一版开发协议的首选 R 初始化；
2. 只读合并两轮既有 G12 实验的 trajectory-derived convergence 结果，审计旧 strict practical convergence 的语义。

本目录不修改任何历史冻结协议、manifest 或结果，不修改 multiFSYNC 模型、先验、CAVI、ELBO 或公共接口默认值，不生成新数据，不启动新 fit，也不运行 continuation。

主要文件：

- `PROTOCOL_R_G12_DEVELOPMENT_ROUTE_V1_20260824.md`
- `R_G12_DEVELOPMENT_ROUTE_V1_20260824.csv`
- `PROTOCOL_CONVERGENCE_SEMANTICS_OFFLINE_AUDIT_V1_20260824.md`
- `audit_existing_convergence_trajectories_20260824_v1.R`
- `convergence_semantics_audit_20260824_v1/`
- `REPORT_G12_CONVERGENCE_SEMANTICS_OFFLINE_AUDIT_20260824_V1.md`
