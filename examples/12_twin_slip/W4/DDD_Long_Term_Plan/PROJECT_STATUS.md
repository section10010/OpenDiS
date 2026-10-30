# 项目状态

更新：2026-10-30；当前 Wave/Gate：W4-G1 准入；W4-G0 PASS（限定有效FCC/I0数值与第二壁机制），G1整体 NOT_RUN。

- 分支 `codex/w4-g0-robustness`，起点 `5b3c900`；不合并 main。
- W1/W2/W3 原始文件保留。规划包 baseline_w3 的19份 Python 与实际 W3 逐字节相同。
- 已现场确认 SSH VM 可达，实际仓库 `/home/sakura/opendis/src`、构建 `/home/sakura/opendis/build`、Python `/opt/miniforge3/envs/anvil/bin/python`。
- W1 原生固定节点补丁仍为14行增、3行删；W3旧副本及W4修改后副本均重新完成A0：10项响应、3项pin、12项原生单元和2项FR回归PASS。
- 24个原始诊断、48个有限细化/消融及8个独立冻结的5 nm修复，共80例真实原生结束，全部数值PASS。
- 新增三阶段速度、反力、接受步线长-时间权重、逐阶段面积与接触账本；真实正例/非法接触反例/CENSORED负控已通过。
- VM参考测试46项PASS，W4配置/谱系/诊断测试8项PASS。Windows缺pytest，参考测试由VM实际执行。
- 20→10 nm尾段2.19%和3.25%的历史UNRESOLVED保留。最后两级主量最大差1.2435%、事件0.7188%，原2%/5%门槛未调整。
- 反向加载先释放、后真实湮灭，HI/HII分别有918/845个非空释放状态，末态塑性闭合通过；首次空数组审查错误已保存、修复，不改原始轨迹。
- 最终审查 analysis_r5：G0 PASS、G1准入PASS；限速不敏感时序与VCoNi仍NOT_QUALIFIED。threader微小非加性UNRESOLVED；loop仅在特定无上限离散路径分辨数值符号，不是材料普适结论。
- 已启动G1真实初始力场小基准：M-only、T-only、M+T共用宏观应力与同一弹性求解器，7项读回/耦合/PK及互反检查PASS；共享接触演化和身份继承尚NOT_RUN，步骤>0显式拒绝。
- VCoNi 材料资格 NOT_QUALIFIED；W4-G1整体及后续均 NOT_RUN，整个W4–W9项目没有完成。

下一步见 NEXT_ACTION.md；所有结果保存到 W4 新目录，不覆盖规划包旧 evidence。G0审查以生产结果而非规划中的NOT_RUN历史标签为准。
