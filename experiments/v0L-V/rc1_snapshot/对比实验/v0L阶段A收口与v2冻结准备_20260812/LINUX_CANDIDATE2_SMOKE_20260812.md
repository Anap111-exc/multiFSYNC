# candidate2 Linux 首次启动微型验收

本说明只运行非正式微型 smoke，不生成 v0L-V 的 10 个正式数据，不启动 210 条正式拟合或 72 条正式审计记录。

前提：服务器检出的项目保持当前相对目录结构，项目根目录示例为 `/home/user/multiFSYNC`；R 及 multiFSYNC 的依赖包已经安装。

```sh
cd /home/user/multiFSYNC
mkdir -p _r_test_lib
R CMD INSTALL --library="$PWD/_r_test_lib" Rcode/multiFSYNC
smoke_root="$(mktemp -d /tmp/multifsync-candidate2-smoke.XXXXXX)"
Rscript "对比实验/v0L阶段A收口与v2冻结准备_20260812/run_v0lv_candidate2_e2e_smoke.R" \
  --output-dir="$smoke_root/result"
test -f "$smoke_root/result/SMOKE_COMPLETE.txt"
grep 'V0LV_CANDIDATE2_E2E_SMOKE_COMPLETE' "$smoke_root/result/SMOKE_COMPLETE.txt"
```

验收条件：命令退出码为 0，`SMOKE_STATUS.csv` 全部为 `PASS`，且完成文件同时声明 `formal_experiment=FALSE`、`formal_data_generated=FALSE`、`formal_fits_started=FALSE`、`formal_continuations_started=FALSE`。

每个模型 fit 内部固定 `n_cpus=1`。服务器正式运行时，只能在相互独立的 fit 任务层调度并行；不得在单个 fit 内再次开多核。

`c2_sha256()` 在 Linux 上优先使用已安装的 R 包 `digest`，否则使用系统 `sha256sum`；Windows 仍可回退到 `certutil`。不维护单独的服务器版模型源码。
