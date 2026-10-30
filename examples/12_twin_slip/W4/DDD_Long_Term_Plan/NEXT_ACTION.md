# 下一动作

1. 从已通过的G1初始力场入口继续：`W4/G1/code/shared_force_probe.py`。不得绕过steps>0的显式保护。
2. 实现单共享CTB的每段host/来源和可释放单侧接触；M/T同一宏观载荷、同一原生弹性网络。
3. 先验证新节点材料侧继承、tag重用反例、非法跨材料碰撞/拓扑的拒绝及账本一致性，再做双侧真实短演化。
4. M-only/T-only/MT消融必须比较相同子网络，不能把两个独立通道动画拼接当共享耦合。针对0.04%量级初始跨侧力明确更严格的误差/因果标准。
5. G1通过后才开发G2周期层片与完整事务回滚；不得直接启动W6集体矩阵。

G0权威运行根：`/home/sakura/opendis/w4_g0_20261030_001/production`；裁决 analysis_r5。
G1运行根：`/home/sakura/opendis/w4_g1_20261030_001`；初始准入 initial_force_probe_r3。
初始准入复跑须新out：`/opt/miniforge3/envs/anvil/bin/python code/shared_force_probe.py --g0-root <G0根> --repo /home/sakura/opendis/src --out <新目录>`。
