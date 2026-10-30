# 项目状态

更新：2026-10-30；当前 Wave/Gate：W4-G1 演化审查；W4-G0 PASS（限定有效FCC/I0数值与第二壁机制），G1整体 UNRESOLVED（矩阵与复核正在执行）。

- 最新原生根 `/home/sakura/opendis/w4_g1_20261030_003`：仅平移力计算原点的数值修复已通过4项原生准入，重新执行完整26例。002旧版互反残差1.68e-6 FAIL保留，原10^-8阈值不变；物理图和材料参数不变。
- G0新A0和四轨迹按既有门禁PASS；额外10^-10完整字段复现FAIL，独立保留，不冒称严格逐字段资格。
- 当前分支 `codex/w4-g1-shared-ctb`，起点为G0提交 `97f0494`；生产入口 `W4/G1/code/run_shared.py`。
- G1运行根 `/home/sakura/opendis/w4_g1_20261030_002`。6项身份几何回归与7项原生准入通过，已真实验证分裂/粗化来源继承、反向单节点释放和跨材料碰撞拒绝。
- 冻结26例演化/消融矩阵；首批14例完整PASS，2.5 nm MT在180秒/0.01821 ns处CENSORED，保留。实测成本后的唯一资源扩展为750秒/例、3600秒总上限，12个新ID从原始初态重跑；物理参数和科学阈值不变。
- 完整反向释放、细化因果审查、原生力反事实及新G0复核全部结束前，不把G1整体标PASS。以下为G0与初始G1历史里程碑。

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
