# OpenDiS W3

本目录是 W3 独立交付，W1/W2 保持原状；工作分支 `codex/w3-source-contact`，没有合并 main。
先阅读 [中文审查报告](REPORT_ZH.md) 或 [带图 HTML](REPORT_ZH.html)。Gate PASS 只针对报告明确的单源和 I0 单通道机制基准；VCoNi 材料资格为 `NOT_QUALIFIED`，W4 为 `NOT_RUN`。

## 目录

- `code/`：唯一配置、FCC 几何、真实原生求解器适配、接触、谱系、回归、运行/收据、统计/制图/报告生成器。
- `results/`：便于审查的派生数据、Gate JSON、配对统计、命令索引、统一 PNG/PDF。
- `evidence/`：完整原始证据，包含失败和被更正的探索；不覆盖历史 result。
- `native_evidence_20261029.tar.gz`：完整原始证据归档。没有对所有文件计算 SHA256。

大型原始证据保留在本地交付目录并由归档同步，不加入源码提交；共享本分支时也需要共享该归档。

## 已核实环境

VM 代码 `/home/sakura/opendis/src`，CMake 构建 `/home/sakura/opendis/build`，Python `/opt/miniforge3/envs/anvil/bin/python`。
原运行根目录 `/home/sakura/opendis/w3_source_contact_20261029`。默认 fish，下面使用 `bash` 语法。
Windows `E:/MiniForge/python.exe` 可运行配置/谱系回归、统计及制图；真实动力学使用 VM 的 ExaDiS 构建。

`evidence/provenance/environment.json` 保存原生模块位置、版本及 Git revision；`w1_native_pin.patch` 是原本就存在的 W1 固定节点补丁，本轮没有改动原生源码。
新环境应核对该补丁并重新构建，通过下方 A0；不要只复制另一个环境的 PASS。

## 新运行环境重新完成 A0

```bash
PY=/opt/miniforge3/envs/anvil/bin/python
CODE=/home/sakura/opendis/w3_source_contact_20261029/code
NEWROOT=/home/sakura/opendis/w3_replay_new
"$PY" "$CODE/qualify_a0.py" --repo /home/sakura/opendis/src \
  --build /home/sakura/opendis/build --root "$NEWROOT"
```

这会在新目录保存精确命令、退出码、配置/几何回归、原生响应、W1 pin 回归和上游 12+2 回归。任何失败都不会生成 A0 PASS。
后续 `run.py` 除了要求 A0 PASS，还逐字节核对三份 A0 冻结生产文件；修改这些文件后必须重做 A0。

## 重放一个正式动力学案例

在已核实的同一 VM 构建中，可以直接用交付的冻结 A0 门槛：

```bash
PY=/opt/miniforge3/envs/anvil/bin/python
ROOT=/home/sakura/opendis/w3_source_contact_20261029
"$PY" "$ROOT/code/record.py" --cwd "$ROOT/code" \
  --receipt "$ROOT/replay_HI_both_receipt.json" --log "$ROOT/replay_HI_both.log" --timeout 360 -- \
  "$PY" "$ROOT/code/run.py" --repo /home/sakura/opendis/src \
  --a0 "$ROOT/A0/gate_result.json" --config "$ROOT/configs/B_loop_HI_both_tau300_t100.json" \
  --out "$ROOT/replay_HI_both"
```

输出目录、收据或日志已存在时会拒绝覆盖。跨机器重放请先使用自己的新 A0 门槛。
初始化重启使用 `--initial 实际构型.npz`，不会重新生成 Burgers；动态释放试验依此实现。

## 重建本轮正式案例（新 ROOT）

下列每一行顺序运行，避免多批任务争用 VM。只执行本轮单源与最小 B：

```bash
"$PY" "$NEWROOT/code/campaign.py" "$NEWROOT" A1 reference
"$PY" "$NEWROOT/code/campaign.py" "$NEWROOT" A1 bracket
"$PY" "$NEWROOT/code/campaign.py" "$NEWROOT" A1 geometry
"$PY" "$NEWROOT/code/extra_checks.py" "$NEWROOT" fine_source
"$PY" "$NEWROOT/code/test_contact_native.py" --repo /home/sakura/opendis/src --out "$NEWROOT/B_contact_unit_final"
"$PY" "$NEWROOT/code/campaign.py" "$NEWROOT" B loop_minimal
"$PY" "$NEWROOT/code/campaign.py" "$NEWROOT" B threader_minimal
"$PY" "$NEWROOT/code/extra_checks.py" "$NEWROOT" B_refine
"$PY" "$NEWROOT/code/extra_checks.py" "$NEWROOT" release
"$PY" "$NEWROOT/code/extra_checks.py" "$NEWROOT" radius
"$PY" "$NEWROOT/code/analyze.py" --root "$NEWROOT" --out "$NEWROOT/analysis"
"$PY" "$NEWROOT/code/figures.py" --root "$NEWROOT" --out "$NEWROOT/figures"
```

`campaign.py` 是本 VM 的便捷入口；其他仓库路径可直接调用 `run.py --repo`。
`B loop` / `B threader`（无 `_minimal`）是原 0.25 ns 探索，不是最小验收所必需，也不是四厚度计算。
本次实际81份收据的命令和退出码见 `results/analysis/all_command_receipts.csv`；A0 的详细原始日志在 `evidence/A0/`。

## 回归与离线分析

```bash
cd "$NEWROOT/code"
"$PY" -m unittest test_configuration test_lineage -v
"$PY" test_contact_native.py --repo /home/sakura/opendis/src --out 新的接触测试目录
```

谱系回归使用图结构夹具测试“带结点脱离”和“碎环并入旧环”这两个记账反例，**不作为真实发射证据，也不用于生产初始化**。

本地重新生成统计和报告：

```powershell
& E:/MiniForge/python.exe code/analyze.py --root evidence --out results/analysis
& E:/MiniForge/python.exe code/figures.py --root evidence --analysis results/analysis --out results/figures
& E:/MiniForge/python.exe code/build_report.py .
```

制图和报告生成器读取指定的最新分析结果及原始网络；不拟合或硬编码厚度趋势。

## 接口边界

当前 I0 为非周期、单凸通道的不可透射接触；全网必须在同一通道内。不能直接用于横跨多个 M/T 片层的周期 RVE。
外部晶体侧、透射/吸收反应、残余 Burgers 和部分位错均未实现。局部交付为线密度/尾段/全局扫掠与塑性闭合，完整空间塑性张量场未完成。
原生 Euler 子步、碰撞、拓扑和 Remesh 后均检查非法穿透；遇到违反约束的状态立即 FAIL，保留失败构型。
NPZ 的 `force_N` 是原生弹性、自力及外应力产生的力，不含独立输出的接触反力；`velocity_m_s` 是原生GLIDE经I0约束和速度上限处理后的实际积分速度。
