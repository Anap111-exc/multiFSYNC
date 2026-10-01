# Stage 6K

这是针对困难数据种子 `84611001` 的一次小型、配对、开发性因果审计。
它只新增 4 个 `G12 + pre-score 1 + Jaoua annealing + gate 10 + fixed 400`
端点，并只读复用 Stage 6J 的 4 个普通 G12 对照。

本目录不改变模型、先验、CAVI、ELBO、公开接口或默认 R 路径；也不授权新数据、
continuation、正式 Monte Carlo 或额外门控长度。完整边界与停止规则见协议文件。
