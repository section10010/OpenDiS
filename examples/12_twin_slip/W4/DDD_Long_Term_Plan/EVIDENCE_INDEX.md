# 证据索引

## 2026-10-30 G1最终裁决（当前权威入口）

- 最终中文报告：`../G1/review/G1_20261030_003/final_review/W4_G1_FINAL_REPORT_ZH.html`；机器裁决`gate_review.json`和`spatial_convergence.json`、`tail_convergence.csv`同目录。
- G1整体UNRESOLVED，尾段2%验收FAIL：M 3.0308416%、T 3.0308825%。1.25 nm MT为NOT_RUN/null；控制器按原预算预留策略CENSORED。
- 实际VM根 `/home/sakura/opendis/w4_g1_20261030_003`。26主例、3原点反控、2细化实际完成；本次只再审查原始输出，没有新增动力学。
- 003完整原始包：`../G1/runs/G1_20261030_003/G1_native_evidence_20261030_003.tar`，**2,603,397,120 bytes，2460文件**。轻量ZIP为20,115,719 bytes、1875文件；仅省略密集accepted_velocity块，完整包保留。
- 002完整失败包：`../G1/runs/G1_20261030_002/G1_native_evidence_20261030_002.tar`，**1,713,602,560 bytes，1897文件**；轻量ZIP为15,540,597 bytes、1406文件。旧力互反FAIL、180秒CENSORED和原报告都保留。
- 本地native视图按ZIP CRC及归档库存尺寸核对，未做大型SHA256扫描。完整归档不是只含构型截图；包含接受步三速度/原生力/反力/dt、配置、命令、源码快照、事件和初末态。
- 本次实际校验收据：`../G1/review/G1_20261030_003/ARCHIVE_VERIFICATION.json`；小型命令/日志/构型逐字节副本在同目录evidence，交付收据为DELIVERY_RECEIPT.json。003 TAR成员核对确认M/T细化密集块存在，未运行的MT结果确实不存在。
- 细化控制器原始命令和退出码：003/spatial_extension_receipt.json（exit 1、wall 1641.98 s）；实际原生M/T均PASS，MT未启动。最终resource_review.json区分这三种事实。
- 003/analysis_r1保留6.34%尾段UNRESOLVED；最终final_review在纳入已有1.25 nm M/T后仍UNRESOLVED，不覆盖首份审查。
- 原生力修复：003/centered_force_admission、force_audit、origin_sensitivity_result.json；002/force_audit、force_invariance为旧失败与定位证据。
- 来源/拓扑：003/admission_r1、admission_r3、identity_equivalence；实际跨材料碰撞76→77后拒绝，first_invalid/restored_legal NPZ在原目录。复杂反应未获资格。
- G0复核：002/g0_regression及003对应逐字节小型副本，原科学门禁PASS；strict_full_field_replay.json额外10^-10标准FAIL保留。密集原始G0复算帧在002完整包。
- 图版：G1/review/G1_20261030_003/reference_figures，均为真实2.5 nm参考的PNG/PDF/SVG；标明最终尾段资格尚未通过，不冒充缺失的1.25 nm MT。
- 生产源码：`../G1/code`；原点修复规格`shared_evolution_centered_v2.json`；本轮最终只读审查器`final_review_shared.py`不调用ExaDiS。
- G2准入BLOCKED，W4-G2与W6集体计算NOT_RUN；整个W4–W9没有完成。

## G0与初始G1历史证据

- 规划准则：DDD_Full_Plan.md、CODEX_NEXT_TASK.md、docs/、configs/gate_criteria.json。
- 真实W3生产与历史：`../../W3/code`、`../../W3/results`；规划副本 baseline_w3 已逐字节核对。
- 新生产代码：`../code`；不改 baseline_w3。新模块 diagnostics.py、g0_campaign.py、g0_analyze.py、g0_extend.py、g0_preflight.py。
- 本轮Windows：`../runs/G0_20261030_001`；原生 VM 新根 `/home/sakura/opendis/w4_g0_20261030_001`，数据在 production/cases，日志在 production/logs。
- 新A0：VM baseline_a0/A0 与 production/A0；新原生接触/反力验收在 production/contact_admission、diagnostic_admission。
- 真实正/反/资源负控：production/preflight_controls.json 和对应构型；参考46项、W4单元8项日志位于production根。
- 判据与原始输入：production/frozen_criteria.json、frozen_inputs；72例清单位于campaign_manifest.json。没有改写原24个输入。
- 审查：production/analysis_r1–r4均保留，以analysis_r5综合裁决为准。Windows在../runs/G0_20261030_001/analysis_r5。
- 中文G0报告和三组真实诊断图：../runs/G0_20261030_001/review_r5/W4_G0_REPORT_ZH.html。
- 随分支保存的实际命令、退出码和小型原生日志：../runs/G0_20261030_001/qualification_evidence/；包括A0、80例命令收据、接触/三速度准入、反向释放及首次审查错误日志。它们从原输出逐字节复制，不替代完整轨迹。
- 完整原始TAR：../runs/G0_20261030_001/G0_native_evidence_20261030_001.tar（3,297,966,080 bytes）；小型native视图已按库存和ZIP CRC核对，密集速度块在完整TAR/VM。
- G1初始真实力场：../G1/initial_force_probe_r3/result.json，三个原生初始网络及力/速度NPZ、唯一规格、宏观张量、24系RSS、局部晶体学读回、原生PK负控和命令收据。R2→R3网络/力回归一致。
- 状态入口：PROJECT_STATUS.md、GATE_TRACKER.md、SCIENTIFIC_BLOCKERS.md、NEXT_ACTION.md。

没有生成新论文图；G0数值/因果PASS、G1初始力场PASS和完整共享接触/材料资格分别裁决。整个长期项目未完成。
