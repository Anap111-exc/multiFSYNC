# Stage 6H reachability causal audit

这是 Stage 6G 真值解封后的最小开发诊断。它复用同一 minimal 数据及 8 个既有
留出端点，新增 16 个留出端点和 12 个完整数据端点。科学方法不变；唯一操纵是
起点数（4 对 12）和训练时间点（5--8 对 6--9）。

服务器入口为 `launch_server_stage6h_20260913_v1.sh`。拟合完成后先冻结三个条件内
的普通 ELBO winner，再执行已授权科学评价。输出写到仓库外，Stage 6G 只读。
