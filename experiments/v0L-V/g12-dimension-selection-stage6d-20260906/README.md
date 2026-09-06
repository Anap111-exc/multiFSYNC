# G12 Stage 6D

本目录是 G12 过指定维度真值盲选择的独立确认实验。它不修改 multiFSYNC 模型、先验、CAVI、ELBO 或 G12 路径，也不是正式论文 Monte Carlo。

冻结设计：3 个全新数据种子、5 个维度配置、每配置 12 个成对随机起点，共 180 条固定 400-sweep 拟合。每条 fit 内 `n_cpus=1`；8 vCPU ECS 上默认 7 个外层 worker。

主要入口：

```bash
Rscript run_stage6d_20260906_v1.R --action=check
Rscript run_stage6d_20260906_v1.R --action=prepare --output-root=/root/g12-stage6d-20260906-v1
Rscript run_stage6d_20260906_v1.R --action=run --output-root=/root/g12-stage6d-20260906-v1
```

`prepare` 生成并密封三个数据，同时建立确定性中间时点留出；`run` 运行 180 条 fit，先按配置内普通 ELBO 冻结 15 个 winner，再进行两个独立维度轴的留出 one-SE 选择。完成后真值仍保持密封，必须等待用户另行授权评价。

服务器脚本 `launch_server_20260906_v1.sh` 会在后台顺序执行 prepare 和 run；`status_server_20260906_v1.sh` 只读报告进度。
