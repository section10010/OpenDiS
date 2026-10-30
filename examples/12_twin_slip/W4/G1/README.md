# G1 共享 CTB 演化入口

当前生产入口 `code/run_shared.py`，规格 `configs/shared_evolution_centered_v2.json`。当前真实运行根 `/home/sakura/opendis/w4_g1_20261030_003`；002为保留的旧版诊断。

002变形末态的原生跨侧力互反残差约1.68e-6，超过原10^-8门槛。003仅将力内核输入原点统一平移到网络附近，实际位错/CTB坐标、材料、载荷、迁移率不变，不投影或配平力；旧失败构型修复后约1.2e-11，4项原生准入通过。完整矩阵和独立原点反控仍在执行，当前整体UNRESOLVED。

材料来源由Burgers、滑移法向及源平面唯一恢复，tag不充当材料身份。每个原生阶段检查并拒绝非法提议。共享CTB在单个弹性网络中对两侧施加相反单侧法向约束，保留切向运动和反向释放。当前仅支持可唯一辨识的仿射单滑移普通环，未知结点反应显式拒绝，不冒充通用多滑移反应。

6项身份几何和7项原生准入PASS。26例动态矩阵与G0复核尚在执行，当前总体UNRESOLVED；首个细网格MT资源CENSORED保留，`resource_extension_r2.json`冻结唯一有界扩展，所有重跑使用新ID，不拼接轨迹。

重放先执行 `qualify_shared.py --g0-root ... --repo ... --spec ... --out admission_r1`，再用 `campaign_shared.py` 在新根冻结并执行矩阵。原生输出目录不可覆盖；全局G0准入与源码字节绑定必须通过。审查工具为 `force_audit.py`、`regress_g0.py`、`analyze_shared.py` 和 `report_shared.py`。

## 历史初始力场入口

已完成真实原生初始共同力场准入 `initial_force_probe_r3`，不是完整G1演化。唯一规格在 configs/initial_shared_probe.json，包含共同应力、迁移率、两侧位置与共享界面。

- 两个FCC取向满足共同[111]的60°关系；同一CTB几何平面 z=2000 nm。
- M-only/T-only/MT使用同一宏观应力张量，同材料参数，单个MT网络在同一ExaDiS全对力求解器中计算。
- 两个单侧子网络与MT对应节点、Burgers和滑移面一致；实际RSS逐段读回，统一σ的原生PK差分与解析单位契约相符。
- 来自另一侧的初始力非零：M/T约5.8724e-12/6.0024e-12 N范数，相对单侧总力0.04457%/0.04556%。跨侧互反残差3.17e-12，PK张量相对误差5.18e-15。

这证明原生网络确实耦合，**没有证明共享CTB接触演化、透射或材料强化**。两侧距CTB尚有间隙，此处只计算初始力。M-only/T-only为消融，不是同密度材料参考。

尚待实现：host/来源继承、tag重用处理、共享接触可释放自由度、跨域碰撞/拓扑一致性、反力事件账本及真实短演化。`--steps` 大于0明确抛出NotImplementedError，不退回两个W3通道。

G1整体状态 `NOT_RUN`。下一动作见规划入口 NEXT_ACTION.md。

```bash
/opt/miniforge3/envs/anvil/bin/python code/shared_force_probe.py \
 --g0-root /home/sakura/opendis/w4_g0_20261030_001/production \
 --repo /home/sakura/opendis/src --spec configs/initial_shared_probe.json --out 新的原生输出目录
```

初次probe元数据曾列出完整设计人口，R2改为仅列活动人口并加入逐段RSS核对；R3移到唯一规格并核对每侧原生滑移面/方向在对应局部晶体中的FCC资格。R2→R3三组实际网络与力已回归一致；旧构型/日志保留，不覆盖。
