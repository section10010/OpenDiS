# 项目状态

更新：2026-10-30。当前 Wave/Gate：W4-G1最终审查；整体 **UNRESOLVED**，尾段2%验收 **FAIL**。本轮报告和归档收尾，不启动新动力学。

- 分支 `codex/w4-g1-shared-ctb`；源码实现里程碑 `4048896`，G0基线 `97f0494`；未合并main。
- 当前原生根 `/home/sakura/opendis/w4_g1_20261030_003`，Python `/opt/miniforge3/envs/anvil/bin/python`，原生仓库 `/home/sakura/opendis/src`、构建 `/home/sakura/opendis/build`。
- 003主矩阵26例完成：M-only/T-only/M+T前向、反向释放、无CTB、原生操作消融；另3例原点敏感性控制完成。
- 1.25 nm细化M与T均数值PASS，分别6954/6938步达到0.025 ns。尾段相对差为3.0308416%和3.0308825%，超过原2%，没有放宽门槛。
- 1.25 nm MT未启动（NOT_RUN）。控制器1641.98 s后剩余约1358.02 s，不足预留完整1500 s单例额度，因此按原3000 s预算策略CENSORED。没有重跑、增额或拼接轨迹；后台原生进程已结束。
- 最终裁决 `G1/runs/G1_20261030_003/final_review/gate_review.json`；最终中文报告同目录 `W4_G1_FINAL_REPORT_ZH.html`，从实际末态NPZ独立复算尾段、守恒和塑性闭合，未运行MT的数据保持null/NOT_RUN。
- 共享单CTB、仿射单滑移来源继承、合法切向运动/释放、阶段检查与非法跨材料反应拒绝已实现并通过限定原生准入。普通闭合环不称为FR源。
- 002旧版力互反残差约1.68e-6 FAIL保留。003仅平移原生力计算坐标，物理图/材料/应力/界面不变，不配平力；4项修复准入及已保存2.5 nm演化构型力审查PASS，原10^-8标准未变。
- 7项身份/能力回归、7项原生准入通过；整数标签优化与610个实际图/133747节点观测逐项等价，新增原生准入r3通过。
- G0生产源码不变；新A0、固定节点/FR、接触/三速度与4条代表性轨迹按原科学门禁PASS。额外10^-10完整字段复现FAIL单独保留。G0原80例R5、W1/W2/W3历史均不覆盖。
- G0限速独立时序、VCoNi标定仍NOT_QUALIFIED；1.25 nm独立时间细化及MT轨迹力审查NOT_RUN。
- G2准入BLOCKED，G2及W5–W9相关未运行项维持NOT_RUN/NOT_QUALIFIED；没有启动四厚度集体计算，整个长期项目未完成。

可审查副本在 `G1/review/G1_20261030_003/`；完整原始归档与小型native视图见EVIDENCE_INDEX。下一恢复动作见NEXT_ACTION；当前不将源码实现或数值完成等同于物理资格。
