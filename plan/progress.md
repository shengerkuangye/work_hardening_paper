# Progress

## 2026-08-24

- Stage: S3→S4（施密特因子与晶粒尺寸证据核对及正文修订）。
- 当前任务包：`plan/task-packets/2026-08-24-schmid-grain-size-manuscript-revision.md`。
- 已确认：15° HAGB 晶粒尺寸未随冷变形量连续降低；2°取向域在中高变形量下变小，但其绝对值受最小像素阈值影响。
- 已重新生成：AD 单轴拉伸几何假设下的六状态 Schmid 因子汇总，临时派生结果位于 `.codex_tmp/manuscript_schmid/06_axial_propensity/`。
- 已修订：`manuscript_draft copy.md` 的方法、EBSD 结果、机制讨论、中英文摘要、关键词和结论；新增表 2 的双尺度晶粒/取向域指标及表 3 的四滑移族 Schmid 因子。
- Review gates: 本次任务的规格符合性审查与论证质量审查均已通过，记录于 `plan/review/2026-08-24-schmid-grain-size-*.md`。

### Capability-use audit

- Required skills: using-research-writing；paper-orchestration；experiment-results-planning；writing-core；verification。
- Skills actually used: using-research-writing；paper-orchestration；experiment-results-planning；writing-core；verification。
- Inputs consumed: `manuscript_draft copy.md`；`docs/2026-07-23-ebsd-mechanism-and-paper-roadmap.md`；`results/mtex_ebsd_comprehensive/README.md`；六状态原始 CTF；重新生成的 `axial_propensity_summary.csv`。
- Inputs not used and why: 未使用 CRSS、滑移迹线、TEM、XRD 总位错密度或 CPFEM 结果；这些数据当前缺失或不属于本次计算，因此只讨论几何取向倾向，不判定实际滑移活动和强化贡献排序。
- Artifacts produced: 修订稿；任务包；表格 schema 更新；施密特因子临时派生 CSV/PNG；规格符合性审查；论证质量审查。
- Verification run: MATLAB R2025a/MTEX 6.1.1 批处理退出码 0，且运行后无残留 `matlab.exe`、`MATLAB.exe` 或 `MATLABWindow.exe`；24 个 raw Schmid 表值逐项核对通过；正文范围检查通过；`style_check.ps1` 退出码 0。
- Remaining risk: 每个状态仅有一个 EBSD 视场，不能估计棒材间或位置间离散性；未给定 CRSS 和旋锻应力路径，Schmid 因子不能识别实际滑移系；项目级 `research_quality_gate.ps1` 因缺少 `chapters/` 和 evidence map 未通过，整篇稿件尚未达到投稿质量门。

## 2026-08-13

- Stage: S3（Results/Discussion 结构与证据规划）。
- 已核对参考论文 PDF 的章节顺序、主要图表及强化模型。
- 已审查本项目拉伸、EBSD、Schmid 因子和强化分解的现有证据边界。
- 当前产物：六模块对照架构及图表/过渡建议（对话交付）。
- Review gates: 规格符合性审查与论证质量审查均已通过，记录于 `plan/review/reference-comparison-*.md`。

### Capability-use audit

- Required skills: pdf；using-research-writing；paper-orchestration；experiment-results-planning。
- Skills actually used: pdf；using-research-writing；paper-orchestration；experiment-results-planning。
- Inputs consumed: 用户提供的 TA16 PDF；用户截图；拉伸趋势表；现有 EBSD 路线图；三项机制审查结果。
- Inputs not used and why: 未逐一读取全部原始 CTF、拉伸原始曲线和 TEM 数据；本轮目标为模块架构，不是重新计算实验指标，且 TEM 尚待补充。
- Artifacts produced: `plan/project-overview.md`、`plan/outline.md`、`plan/progress.md`、`plan/task-packets/2026-08-13-reference-comparison.md`、`plan/review/method-experiment-traceability.md`、`plan/experiment-protocol.md`、`tables/table-schema.md`、`figures/data-manifest.md`。
- Verification run: 对参考 PDF 16 页进行文本提取并核对 3.1–4.3 小节、图 3–16 和强化计算；与本地拉伸及 EBSD 汇总交叉核对。
- Remaining risk: 本炉成分、总位错密度、TEM 滑移类型、体织构及旋锻多轴载荷路径尚缺；Rp0.2 存在实验室记录与原始曲线复算口径差异。
