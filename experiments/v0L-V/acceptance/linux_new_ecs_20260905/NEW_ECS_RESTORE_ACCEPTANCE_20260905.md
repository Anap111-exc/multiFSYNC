# 新 ECS 系统盘恢复与 Linux 验收（2026-09-05）

## 结论

新实例 `i-bp1egmsairexm9w1rnks` 已由初始公共镜像切换到自定义镜像
`m-bp11iqykawepnl38huto`。恢复后的源码、依赖、历史实验目录和 frozen1
绑定均可用。安全更新后的 candidate2 端到端微型 smoke 为 12/12 PASS，
G12 公共接口 smoke 通过。本次未启动正式数据、正式 fit、正式 continuation
或正式评价。

## 实例与系统盘

- 规格：`ecs.u2a-c1m2.2xlarge`，8 vCPU，控制台标称 16 GiB。
- 系统：Alibaba Cloud Linux 3.2104 U13.3，x86_64。
- R：4.5.3；Git：2.43.7。
- 根设备：`/dev/nvme0n1p3`，ext4；系统盘约 50 GB。
- 验收时根分区约 49 GB，已用 20 GB，可用 28 GB，使用率 42%。
- Linux 显示可用内存总量约 14 GiB；无 swap。
- 恢复后的 Stage 6A、6B、6C 服务器目录均存在。

## 源码与 Git 完整性

- 恢复归档中 candidate2 的 50 个运行绑定文件全部存在，SHA256 为 50/50。
- 因服务器直连 GitHub 的 Git HTTPS 仍超时，本机从已核验仓库制作完整 Git
  bundle，经 Workbench 上传；bundle SHA256 为
  `b94bba7097d76de94204a9c146f2234a622dff1481357971df7067a462c1b750`。
- 服务器独立检出目录：`/root/multiFSYNC-frozen1-git`。
- 分支：`release/v0L-V-rc1`。
- HEAD：`aee98b79a80a6535a0ad7a7c5187b29db2f66176`。
- 已存在的标签 `v0.3.0-v0L-V-frozen1` 解引用到同一 commit。
- 工作树干净；`git fsck --full --strict` 通过；远程 URL 指向正式 GitHub 仓库。
- frozen1 自带验证器结果：静态检查 PASS、运行源码 SHA256 50/50、字节数
  50/50、manifest 10/210/72、运行时检查 PASS、正式执行未授权且未开始。

## 安全更新后 smoke

- 安装 21 个安全更新包；更新后待处理安全条目为 0。
- `needs-restarting -r` 判定无需重启。
- candidate2 完整非正式微型链：12/12 PASS，退出码 0，墙钟 4.17 秒，
  最大 RSS 194672 KiB；stderr 为空。
- G12 公共接口：PASS，退出码 0，墙钟 0.22 秒，最大 RSS 86064 KiB；
  stderr 为空。
- candidate2 COMPLETE 记录明确为：
  `formal_data_generated=FALSE`、`formal_fits_started=FALSE`、
  `formal_continuations_started=FALSE`。

## 短并发预检

每个进程均固定 BLAS/OpenMP 线程为 1；并发只发生在独立任务层。

- 3 个独立 candidate2 e2e 微型任务：3/3 PASS，墙钟 4.450 秒；各任务
  最大 RSS 的和为 590312 KiB，单任务最大值 197400 KiB。
- 6 个独立 candidate2 e2e 微型任务：6/6 PASS，墙钟 5.916 秒；各任务
  最大 RSS 的和为 1177952 KiB，单任务最大值 197520 KiB。

这里的 RSS 和是各进程各自峰值的保守上界，不等于同一时刻实测总峰值。该
测试确认 3/6 并发调度和线程隔离正常，但微型 smoke 不能单独替代全尺寸拟合
的内存与运行时间证据。

## 边界

- 未覆盖旧的 `/root/multiFSYNC-rc1` 归档目录；精确 Git 检出使用新目录。
- 未修改模型、先验、CAVI、初始化、退火、停止或评价公式。
- 未生成正式数据；未启动正式 fit、continuation 或评价。
