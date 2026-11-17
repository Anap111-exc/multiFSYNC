# multiFSYNC 原CAVI实验冻结备份 · 2026-11-17

本分支保存 `cavi_true222_fit444_m4_20261117_v1` 的隔离执行源码、注册配置、必要测试、来源证明和哈希。它是独立快照分支，不替换原项目分支，也不改变已在ECS运行的实验。

## 本次实际方法

- 原CAVI；真实共享／两个研究特异因子数为(2,2,2)，拟合容量(4,4,4)。每个因子真实FPCA秩2，拟合容量4。
- 2个研究，每研究60人，500变量；10个新总体，各12起点，共120条正式拟合。
- random＋gram_unit_energy＋一次预得分初始化；99步退火＋400步普通更新。
- 同一数据集内按最大有效完整普通ELBO选择；PPI≥0.5计数，选择时不读取真值。
- 源码兼容MGP接口，但本次 `factor_shrinkage=NULL`／`mode="none"`，未启用MGP。

## 文件与执行版本

`snapshot/` 保留原相对路径和精确文件字节。原108文件清单及两份总清单完整保存；另外保留修订、相关测试和小型启动验收证据。

1. `snapshot/CODE_SHA256SUMS.txt` 与 `PACKAGE_MANIFEST.json`：最初上传包，原样保留。
2. `snapshot/amendments/probe_call_v2/`：独立两步资源探针的显式NULL调用修复。第一次调用在参数检查阶段失败，未执行CAVI；正式120条拟合设置不变。
3. `snapshot/amendments/ops_v3/`：实际运营控制入口，存储预算6.75GiB＋3GiB余量，以及控制器内部的首端点观察线程；原数值源码不变。
4. `snapshot/amendments/resource_observer_v1/`：观察器及测试依赖。独立进程版仅保留作测试依赖；实际运行线程版。
5. `snapshot/evidence/RUNTIME_ACCEPTANCE.json` 是云端 `RUNTIME_ACCEPTANCE_V3.json` 的同字节本地副本；`REGISTERED_LAUNCH_VERIFIED.json` 仅表示启动时核验，不是当前进度或科学结果。

实际云端入口为 `amendments/ops_v3/launch_registered_pipeline_v3.py`，数值入口仍为原执行根目录中的 `scripts/run_cloud_fit.R`。配置、路径及调用关系详见 `snapshot/LAUNCH_REPORT.md`，不能把早期探针脚本当作正式拟合入口。

## 核验

运行 `python verify_backup.py`。它只读取文件并校验总清单、原执行清单及修订哈希，不启动R拟合、不连接ECS、不调用云API。

Git行尾转换已关闭。SHA256反映原文件字节，不能在未更新版本记录时统一换行或格式化。备份提交之后还须逐文件核对Git对象及远端取回副本，验证记录留在本地备份回执。

## 恢复边界

这是一份源码和配置快照，不是单仓库便携安装包。云端流程依赖已冻结的Linux multiFSYNC安装库、Splinets及其他R包；依赖位置和文件哈希见 `RUNTIME_LIBRARY_BINDING.json`、`scripts/validation/HOST_INSTALLED_PACKAGE_HASHES.csv`。快照不包含安装库二进制、凭据、实际数据、真实参数、拟合端点或历史归档，完整实验恢复需另行取回并核对相应数据和Linux环境。

`generation/test_local_contract.R` 与 `generation/prepare_registration.py` 引用了历史E盘来源，是历史验收／注册溯源记录；不能声称在任意新电脑上直接运行。正式实验使用已经冻结的CSV，无须重新注册或更换种子。Linux冻结数据绑定的基、结点和变换不能用Windows重新生成版本替代。

## 后续分支使用

本备份分支固定保留，后续复现锁定本次提交。后续算法或正式模拟框架开发，应从经过核验的提交另建开发分支，完成验证后为每轮实验记录新的冻结提交。不能把修改后的开发分支继续称为本次实验原版本。

本分支为独立根提交的备份树；其他历史版本仍保留在仓库的既有分支中。本次备份未切换、提交或清理D盘工作树，未干预运行中的ECS任务。
