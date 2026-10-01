# Stage 7L 执行记录（启动阶段）

2026-10-01 12:38（北京时间）已通过 Workbench 启动。
当前记录仅确认启动与核验通过；尚未获得科学结果。

## 版本与核验

- 本地 Git 分支 dev/v0lv-init-geometry-20260818，HEAD 5ec35267def4ccbc1a71c6811a9b06acccb19e9c。用户原有 dirty 修改保持。
- 用户给定生成器、pilot、Stage7C runner、UAS helper 的四项 SHA256 均匹配。
- 沿用 Stage7I 安装包；源码归档及安装包 DESCRIPTION/NAMESPACE/rdb/rdx 全部匹配。DGP 的 Splinets 使用已有 /root/multiFSYNC-jaoua-paired-lib-20260909，安装版本和哈希另行登记；Stage7I 库优先，未覆盖历史包。
- 引用既有 Stage7J 全449步复现：19307项状态检查、400项普通ELBO，差值均0。未重复长重放。
- 新适配本地微型检查只续跑既有toy的51、52，两步原始状态逐字段通过，逐快照落盘通过。
- ECS源码校验、R解析、评价器绑定及 RAM StopCharging dry-run 通过。R 4.5.3，8核；首批8个 R 进程各约100%单核。
- 新5个C1数据已生成、封存及QC通过；首批U拟合已开始。75个数据/真函数/载荷/均值/beta/fit种子在生成前登记，无试跑筛选。

## 调度和预算

50U（每数据10条）后冻结各U winner，再构造固定S并续跑最多5条，最后统一冻结和评价。8外层worker，单fit及BLAS/OpenMP均单线程。没有A拟合或密集逐sweep状态。

拟合预计10–12小时，即10月1日22:38至10月2日00:38左右；排期预留至10月2日02:38–04:38。技术上限：单fit3小时，整批拟合18小时，自首批U启动计。科学评价及报告由同一自动流水线在盲层冻结后执行。

完成后无其他R任务则用既有RAM profile multifsync-self-stop请求StopInstance(StoppedMode=StopCharging, ForceStop=false)，先dry-run再实调用，并保留响应。系统运行失败也保存状态并尝试同一节省停机；存在其他R任务时跳过停机。

## ECS 位置

- 实例：i-bp1egmsairexm9w1rnks，cn-hangzhou；启动时公网IP 47.98.87.163。
- 执行源码：/root/stage7l-exec-20261001-v1/
- 主systemd unit：stage7l-prospective-20261001-v1.service
- 调度日志：/root/multifsync-results/stage7l_supervisor_20261001_v1/RUN.log
- 本批结果：/root/multifsync-results/stage7l_prospective_20261001_v1/
- 盲层冻结：上述结果目录 blind_freeze/，必须有 BLIND_ENDPOINTS_FROZEN.txt。
- 科学小表：scientific_evaluation/PAIRED_DATASET_SCIENCE.csv、ALL_ENDPOINT_SCIENCE.csv、ALL_ROLE_SCIENCE.csv、PAIRED_ROLE_DAMAGE.csv。
- 报告：REPORT_STAGE7L_PROSPECTIVE_20261001_v1.md。
- 全部哈希：COMPLETE_FILE_HASHES.csv；完成标记 EXPERIMENT_COMPLETE.txt。
- 最终小表归档：stage7l_small_readout.tar.gz，不含原始RDS。
- 节省停机记录：/root/multifsync-results/stage7l_autostop_20261001_v1/。

## 本地保留与后续读取

本地只有协议、清单、源码、微型检查及 preflight_ecs_v1 小表。正式fit、snapshot50、S快照及sealed truth留ECS，不下载大状态。

后续先查完成/失败标记，再下载最终小表归档并核对哈希；ECS已节省停机时须在实例恢复运行后读取。若冻结/接口/源码门槛失败，不重新生成数据或补拟合，报告保留证据。
