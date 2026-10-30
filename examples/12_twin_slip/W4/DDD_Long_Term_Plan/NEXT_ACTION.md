# 下一动作

本轮已完成最终审查、报告与归档收尾；没有后台原生计算。当前G1整体UNRESOLVED，G2准入BLOCKED。不从初始力场检查或26例主矩阵重新开始。

1. 先读 `G1/review/G1_20261030_003/final_review/W4_G1_FINAL_REPORT_ZH.html`、`gate_review.json`及`spatial_convergence.json`。M/T尾段约3.03%仍超原2%；MT细化NOT_RUN，控制器CENSORED证据保留。
2. 下一科研修复应在固定1.25 nm空间网格下设计独立时间步/节点位移敏感性，材料、宏观σ、初始图及I0规则保持不变；先冻结预算和原2%门槛，不以继续加密直到碰巧通过为目标。
3. 增加阶段前后接触端点/尾段审计，区分真实运动与Remesh局部效应；保留原尾段验收，辅助带宽统计不得替换失败观测量。
4. 现有数据先使用；不重复已完成的原生矩阵，不拼接CENSORED轨迹，不将缺测填零。若以后开展新计算，须使用新的case ID及明确输入/资源约束。
5. G1严格资格通过前不启动G2周期层片或W6四厚度集体矩阵。

实际VM：`/home/sakura/opendis/w4_g1_20261030_003`。现有结果再审查入口（只读已有输出，不运行ExaDiS）：

```bash
/opt/miniforge3/envs/anvil/bin/python code/final_review_shared.py --root <003根> --out <新的审查目录>
```

002完整失败数据、003修复版完整数据、W3历史和G0 R5均保留；源码在独立 `codex/w4-g1-shared-ctb`，不合并main。
