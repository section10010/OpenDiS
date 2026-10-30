# W4 执行入口

继续工作先读 `DDD_Long_Term_Plan/PROJECT_STATUS.md`、`NEXT_ACTION.md` 和最新审查报告。
规划原文、W3及旧规划 evidence 保留原状；新生产代码在 `code/`。

## W4-G0 已执行

真实VM根 `/home/sakura/opendis/w4_g0_20261030_001/production`，Python `/opt/miniforge3/envs/anvil/bin/python`，原生仓库 `/home/sakura/opendis/src`、构建 `/home/sakura/opendis/build`。这些路径本轮已现场核实。

- 24原始输入、48有限敏感性与8个预冻结空间修复：80例原生计算，均数值PASS。
- 第一轮20→10 nm的两项尾段未达到2%，保留UNRESOLVED；5/2 nm修复及时间独立减半后，最后两级主量最大差1.2435%、事件0.7188%，未放宽门槛。
- 4000 m/s限速改变部分固定时间面积约2%–7%，微小非加性翻转符号。实际时序与材料预测仍NOT_QUALIFIED。
- 原始/接触后/截断后三速度、反力、真实接受dt和半连接线长权重均保存并复算。
- 原生正例、非法接触反例、资源CENSORED负控、两个实际反向释放和末态塑性闭合通过。
- 新W3及修改后W4分别重做A0；W1 pin补丁14增3删保留，未更改上游C++。

最终报告 `runs/G0_20261030_001/review_r5/W4_G0_REPORT_ZH.html`；最终数值裁决 `analysis_r5/gate_review.json`。
R1/R2及早期R3/R4审查保留；以包含动态释放完整审查的R5为准。
报告中的图是实际G0诊断，不是新论文F1–F10或四厚度结果。

## 重放

新根先运行下列 A0 和原生准入；`inputs/` 提供随分支交付的最小配置来源，原规划包保留原状。先将本目录的路径设为 `W4`。

```bash
PY=/opt/miniforge3/envs/anvil/bin/python
ROOT=/path/to/new/root
W4=/path/to/checkout/examples/12_twin_slip/W4
PLAN="$W4/inputs"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
"$PY" "$W4/code/qualify_a0.py" --repo /home/sakura/opendis/src --build /home/sakura/opendis/build --root "$ROOT"
"$PY" "$ROOT/code/test_contact_native.py" --repo /home/sakura/opendis/src --out "$ROOT/contact_admission"
"$PY" "$ROOT/code/test_native_diagnostics.py" --repo /home/sakura/opendis/src --out "$ROOT/diagnostic_admission"
"$PY" "$ROOT/code/g0_campaign.py" --root "$ROOT" --plan "$PLAN" --repo /home/sakura/opendis/src --prepare-only
"$PY" "$ROOT/code/g0_preflight.py" --root "$ROOT" --repo /home/sakura/opendis/src
"$PY" "$ROOT/code/g0_campaign.py" --root "$ROOT" --plan "$PLAN" --repo /home/sakura/opendis/src
"$PY" "$ROOT/code/g0_analyze.py" --root "$ROOT" --out "$ROOT/review_base" --audit-velocities
# 只有实际主量不收敛才调用；不放宽标准，最多16例/1800秒。
"$PY" "$ROOT/code/g0_extend.py" --root "$ROOT" --review "$ROOT/review_base" --repo /home/sakura/opendis/src
"$PY" "$ROOT/code/g0_release.py" --root "$ROOT" --repo /home/sakura/opendis/src
"$PY" "$ROOT/code/g0_analyze.py" --root "$ROOT" --out "$ROOT/review_final" --audit-velocities
"$PY" "$ROOT/code/g0_report.py" --root "$ROOT" --analysis "$ROOT/review_final" --out "$ROOT/report"
```

`g0_campaign.py` 冻结门槛和全部候选配置，顺序执行而不改原24个输入。每例不可覆盖，遇原生数值FAIL停止依赖计算；资源达到上限保存checkpoint并标CENSORED。精确累计账本重启仍NOT_RUN，不能直接拼接分段当连续轨迹。

`--prior-velocity-audit 旧审查目录` 可复用同一不可覆盖case的已完成三速度审查，但核对状态、终点、步数和面积完全相同并记录来源。它不绕过新case检查。

## 数据与范围

`accepted_velocity` 的chunk含实际左端状态三速度、位置、tag、原生力、反力、乘子、权重与真实接受dt。用于计数的step是其后的接受步；左端time和offsets明确给出。帧文件只是构型显示，不能冒充逐步塑性。

主面积字段为积分线段四边形的累计无符号扫掠；另给出由原生全局塑性独立提取的单滑移有符号等效面积。全局拓扑面积张量逐阶段记录；完整空间拓扑塑性场尚未实现。无弹性势能API时保持NOT_RUN，不把由力导出的功改名独立势能。

完整原始数据原位保留VM，同时交付完整TAR；小型view ZIP/展开目录只省略密集 `accepted_velocity`，其余日志、构型、检查点和代码均在。需要本地 `--audit-velocities` 时先展开完整TAR，或在VM运行。

`code/analyze.py`、`figures.py`、`campaign.py` 等无g0前缀的旧工具是W3兼容回归，不作为新G0入口。未实现周期/I1/材料功能不得通过配置静默开启。
