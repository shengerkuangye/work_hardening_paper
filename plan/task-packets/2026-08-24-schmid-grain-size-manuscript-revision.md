## Task Packet

- Scope: 在 `manuscript_draft copy.md` 的方法、EBSD 结果、机制讨论、摘要和结论中补充施密特因子与晶粒尺寸内容。
- Stage: S3（真实 EBSD 派生结果核对）→ S4（正文修订）→ S5（质量与证据核验）。
- Files to read: `manuscript_draft copy.md`；`docs/2026-07-23-ebsd-mechanism-and-paper-roadmap.md`；`results/mtex_ebsd_comprehensive/README.md`；`.codex_tmp/manuscript_schmid/06_axial_propensity/axial_propensity_summary.csv`。
- Files allowed to edit: `manuscript_draft copy.md`；`plan/progress.md`；`plan/task-packets/2026-08-24-schmid-grain-size-manuscript-revision.md`；`plan/review/2026-08-24-schmid-grain-size-*.md`；`tables/table-schema.md`。
- Required skills: using-research-writing；paper-orchestration；experiment-results-planning；writing-core；verification。
- Evidence/data inputs: 六状态 raw EBSD 的 2°取向域与 15° HAGB 晶粒形貌汇总；在 AD 单轴拉伸假设下四类滑移族的面积加权平均最大绝对施密特因子。
- Required artifacts: 修订后的正文；更新后的表 2（晶粒/取向域尺寸）与表 3（施密特因子）；规格符合性审查；质量审查；capability-use audit。
- Rejection checks: 不将 2°取向域称为物理晶粒；不声称 HAGB 晶粒连续细化；不以施密特因子单独判定实际滑移启动或强化贡献；不将单视场内晶粒/像素当作独立实验重复；不覆盖原始实验数据。
- Validation commands: 核对表内数值与派生 CSV/路线图；检查 `git diff --no-index -- manuscript_draft.md "manuscript_draft copy.md"`；运行学术写作风格检查；扫描禁用词、占位符、表号和关键限定语。
