# Grade 4 工业纯钛旋锻冷变形项目

## 当前状态

本仓库研究 Grade 4 工业纯钛棒材在室温旋锻冷变形过程中的拉伸响应、加工硬化行为和显微组织演化。当前编辑稿为 [manuscript_revision_2026-09-08.md](/home/abcd/repos/work_hardening_paper/manuscript_revision_2026-09-08.md)，按力学性能、晶粒与晶界、ODF/Schmid、TEM 预留、GND 及强化贡献展开，目标为中科院一区、具体期刊待定。本轮已完成可由现有资料支持的全文修订；需人工逐项确认的内容见[模块审核顺序](/home/abcd/repos/work_hardening_paper/docs/2026-09-08-manuscript-module-review.md)及[数据与方法审核记录](/home/abcd/repos/work_hardening_paper/docs/2026-09-08-manuscript-evidence-audit.md)。

[2026-09-06 初稿](/home/abcd/repos/work_hardening_paper/manuscript_framework_draft_2026-09-06.md)作为本轮比较基线保留，[manuscript_draft.md](/home/abcd/repos/work_hardening_paper/manuscript_draft.md)为更早历史稿。旧口径图件、旧强化模型和编写注记见[支持记录](/home/abcd/repos/work_hardening_paper/docs/2026-09-08-manuscript-supporting-records.md)，不作为当前正文的已验证结果。需要 Word 版本时应由当前编辑稿生成。

六个研究状态的名义直径为 7.00、6.5、6.02、5.60、5.25 和 5.00 mm，对应名义截面减缩率 0、13.78%、26.04%、36.00%、43.75% 和 48.98%。旧文件中的 6.48 mm／14.31% 保留为历史标记，原始数据不重新缩放。中文材料通称统一为“工业纯钛”，英文采用 commercially pure titanium；材料牌号、内部标识与历史写法的说明见 [材料命名与内部追溯记录](/home/abcd/repos/work_hardening_paper/docs/material_naming_and_traceability.md)。

## 目录职责

- `data/`：原始拉伸数据、CTF、金相源文件和材料质量证明。原始文件只读使用。
- `tools/`：拉伸绘图、图像拼接及 MTEX/MATLAB 分析脚本。
- `results/`：可由脚本重新生成的分析结果和论文候选图。
- `figures/`：当前正文直接引用的图和对应派生表。
- `references/`：正式文献原文及引用库。
- `extracted_literature_notes/`：文献提取笔记和适用边界。
- `docs/`：分析设计、交接说明和 Linux MTEX 运行规定。
- `plan/`：任务历史、审查记录和尚待补充的数据项。

## 已具备的数据

- 拉伸：六种直径、每组两条原始曲线；6.02 mm 的 Y-2 屈服段保持不可靠标记，当前汇总按人工排除表处理。
- EBSD：六个状态的 raw/denoised CTF 配对，共 12 个 600×600 扫描，步长 0.5 μm；raw 为定量主来源，denoised 仅用于敏感性比较。
- 金相：六个状态的纵截面 100×、200×、500× 图像及源 PPTX。
- 材料证明：`data/material_certificate/gr4b23271_grade4_quality_certificate.png`。其中批次对应关系仍应在投稿前与试样台账复核。
- MTEX：77 个 MATLAB 脚本，覆盖晶粒/取向域、晶界、KAM/GROD/GOS、织构、c 轴分布和 Schmid 因子等分析。
- OIM 拼图：Python 脚本已支持 Windows Arial/Calibri 和 Linux Liberation Sans/Noto Sans/DejaVu Sans；可用 `OIM_FONT_DIR` 指定字体目录。
- Fig. 7 参考式织构图：六状态 PF/IPF 小图及组合大图位于 `results/mtex_fig7_texture_gallery/reference_style/`；状态独立色标版用于观察峰位和形态，原统一 0–7 MRD 版用于跨状态比较。

## 当前主要结论边界

- 抗拉强度随冷变形量整体升高，断后延性整体降低；屈服强度存在非单调变化，后续均匀变形区间和加工硬化曲线尚需复核。
- 15° HAGB 晶粒尺寸没有连续降低；现有结果支持晶粒沿 AD 拉长和晶内低角度取向域细分，不支持把强度变化简单归因于连续物理晶粒细化。
- KAM、GROD、GOS 和 LAGB 是位错相关取向梯度或变形亚结构的代理量，不等于总位错密度。
- Schmid 因子只表征指定载荷下的几何有利度，不能单独判定旋锻过程中的实际滑移系。
- 现有 KAM 在 2°重构取向域内计算，排除跨域邻点；各状态有限 KAM-GND 像素的覆盖比例不同，比较需考虑统计区域差异。
- 每个状态目前只有一个 EBSD 视场，不能把像素或晶粒数量作为独立实验重复。

## 仍需完成

- 核对质量证明与实际试样批次，完善成分、旋锻道次、速度、温升和润滑记录。
- 增加独立加工棒材和多位置 EBSD 重复，重点复核 26.04%状态。
- 如需定量强化分解，补充 XRD/TEM 或其他直接位错证据及不确定度传播。
- 统一正文全部参考文献并按目标期刊格式整理；当前已正式纳入 Xia 等 2026 年 TA16 冷轧研究作为方法与机制比较文献。
- `results/mtex_ebsd_comprehensive/README.md` 描述的是完整输出合同；当前仓库仍以分模块结果为主，不能把该 README 的清单误认为所有综合产物均已生成。

## Linux MTEX 运行

所有 Linux MTEX 批处理必须通过：

```bash
tools/run_mtex_guarded.sh --timeout 7200 "MATLAB expression"
```

完整规定见 [`docs/mtex-linux-workflow.md`](docs/mtex-linux-workflow.md)。不要直接后台运行 `matlab ... &`。
