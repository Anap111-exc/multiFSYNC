# Stage 7C 本地实现验收

日期：2026-09-27  
状态：`LOCAL_ACCEPTANCE_PASS`

## 验收范围

本次只验证 Stage 7C 对已接受 36-fit Route-B 含协变量运行时的最小调度扩展。模型、先验、CAVI、ELBO、G12、预得分、退火、400-sweep endpoint、objective eligibility 与真值盲选择规则均未改变。

## 结果

- runner 与定向测试脚本解析：PASS；
- protocol、data manifest、fit manifest 和用户授权 SHA256 绑定：PASS；
- 运行时 registration：9 个数据臂、108 个唯一 fit、每 seed-arm 12 个起点：PASS；
- 每个逻辑数据内 C0/C1/C2 fit seeds 一一配对：PASS；
- 所有候选执行与真值访问字段仍为 `FALSE`：PASS；
- 本地端到端微型 smoke：3/3 fit COMPLETE 且 objective-eligible；
- 微型 truth-free selection：3/3 臂各冻结 1 个 winner；
- 定向回归测试：22/22 checks PASS；
- warning 被保存在 terminal record 与 `WARNINGS.csv`，没有被丢弃；
- smoke 使用独立 seeds `99719001--99719005` 与 `99719601`，未使用三个正式候选 seed 组；
- smoke 未启动科学评价、未解封真值；
- `git diff --check`：PASS。

微型拟合出现预期的 `Max iterations reached before convergence` warning；smoke 只运行极短迭代，该 warning 被完整保留，不影响接口验收。

## 继承源码核验

本轮没有修改 multiFSYNC 包源码。关键源码 SHA256 与已经运行36-fit先导的绑定一致：

- `R/generate_data_structured.R`: `226131949adf52a28fc5e8e91b765a8cdb98a4e2c3a2ffbfa7ee40fc7bc932dd`；
- `R/update_beta.R`: `7194257f581f48338590c595b10121f0cd2186820513bc6ddfb1d7bc425e2876`；
- 继承的36-fit运行时：`5eaa571692b497f3abbf0c4a3083edaa11d9547d6d7f24d972c9a14e2458729d`。

因此没有重复运行完整 testthat；本次新增逻辑由 registration check、端到端 smoke 和22项定向检查覆盖。

## ECS 启动边界

本地验收满足用户授权中的“通过后启动 ECS”条件。服务器仍须依次通过：源归档 SHA256、关键文件逐项 SHA256、独立包安装、Linux registration check、Linux 微型 smoke、正式数据配对与密封 QC、全尺寸 benchmark。只有全部通过后才可启动108条拟合；真值解封和科学评价仍未授权。
