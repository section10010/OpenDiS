# 证据索引

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
