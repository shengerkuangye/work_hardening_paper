# 全文证据与可复现性审查记录

日期：2026-09-08。用户目标为中科院一区，期刊待定。当前结论是：已具备组织与性能比较的基础，但存在影响机制和贡献结论的实际数据/方法缺口；仅靠文字润色不能视为已达到投稿定稿状态。

审查基线为[2026-09-06 稿](/home/abcd/repos/work_hardening_paper/manuscript_framework_draft_2026-09-06.md)，修改结果见[2026-09-08 稿](/home/abcd/repos/work_hardening_paper/manuscript_revision_2026-09-08.md)。下列等级用于安排科学核查优先级：P0 表示直接影响主要结论或可复现性，P1 表示影响解释精度、代表性或投稿图文质量，P2 表示编辑整理。

## 核心发现与处理

| ID | 等级/模块 | 发现及依据 | 本轮已经处理 | 仍需完成 |
|---|---|---|---|---|
| A01 | P0 / M01 | 厂家、INSTRON 5582、GB/T 228.1-2010、无中间退火、Kroll 和电解抛光等仅继承旧稿；旧编写注记 B1 已说明来源层级 | 从已确认方法叙述中移入显式核对项；保留尺寸、标距和 CTF 采集参数 | 原始工艺、设备、标准版本及制样记录；M 初态历史及试样/材料证明对应 |
| A02 | P0 / M02 | Rp0.2 汇总具有拟合字段，但本轮有界检索未找到生成拟合区间和偏移交点的原始程序 | 保留既定曲线口径，明确汇总核对不等于提取验证，局部回落不直接解释为软化 | 恢复或建立可审查提取流程，与机台记录逐试样比较 |
| A03 | P0 / M02 | 曲线与旧实验室 Rp0.2 在 5.60、5.00 mm 状态分别相差 −55.01、−70.94 MPa；6.02 还存在统计集合差别 | 旧性能图、旧模型与正文分开，全部当前性能使用同一转写表 | 查清差异原因，禁止为改善拟合而切换强度基准 |
| A04 | P0 / M02 | 12 条曲线满足 Engineering Strain=DeformValue/50；字段关系不能确定物理测量通道，部分起点高应力/末端应变异常 | εe,Rm 保持“记录峰值位置”定义，不写为 Ag、A 或直接均匀延伸率 | 核对通道、零点、预载、切换及卸载过程 |
| A05 | P0 / M02 | 旧图 3/4 截止点取最大换算真应力位置、负导数截零；平滑窗口还需核对 | 移至历史图清单，题名和主结果不承诺已完成加工硬化阶段分析 | 确认有效均匀区间、应变定义、平滑/求导及所有有效重复；本轮不重画 |
| A06 | P0 / M03、M06 | KAM 在预先赋予 2°重构 grainId 后计算；MTEX 会排除跨域邻点，并非仅由一阶邻域和 5°阈值决定 | 方法、图注、GND 解释补充域内限制和有效像素口径 | 核对掩膜及选择效应，决定是否需要统一处理后的敏感性分析 |
| A07 | P0 / M08 | 旧强化模型用旧强度、固定 M、文献参数及未完整对应的成分；模型的三项不是晶粒/取向/位错三类独立贡献 | 旧 4.4 数值、表 4、图 15/16 整体转入追溯文件；正文改为综合讨论 | 在统一强度基础上约束参数、取向和密度；给出不确定性或保持非定量结论 |
| A08 | P0 / M09 | Meng 等已经研究纯钛旋锻的组织、织构和强化，三因素联合分析不是已证实研究空白 | 引言增加直接文献对话，收紧为具体组织尺度与亚结构问题 | 扩展最接近研究的系统比较，确认新增认识具有代表性 |
| A09 | P1 / M03、M07 | 纵截面 ECD 不变不意味着各方向边界间距不变；LAGB 长度分数不等于边界长度密度 | 4.1 修正为基于当前二维 ECD 的解释范围，不排除其他晶界几何作用 | 核实形状各向异性、取样方向；需要定量归因时补相应量 |
| A10 | P1 / M04 | 原 Schmid 图调用 pyramidalCA，其定义为一阶锥面〈c+a〉；高 SF 或 b 类型不确定实际滑移面 | 正文及图注明确一阶锥面，CRSS 仍为文献输入 | 核对轴向映射、组分和滑移迹线/TEM；不由两个单晶给出织构百分比 |
| A11 | P0 / M03、M06 | 每状态只有一个 EBSD 视场；有效像素掩膜随状态变化 | 限定为所测视场与给定处理条件，不以像素/晶粒数量作实验重复 | 核对径向位置和独立加工重复；按目标差异设计补充取样 |
| A12 | P1 / M05 | TEM 未开展、MD 仅有方案 | 保留 TEM 方法和结果占位，不写预期观察或已做 MD | 先确定要辨识的局部组态，再按实际结果决定是否需要 MD |
| A13 | P1 / M11 | 旧力学/组织图仍有 6.48/14.31 标记；金相图方向及比例尺、IPF 色键完整性待核 | 图件索引明确现有文件状态；只重编号正文，不改原图 | 投稿前按核实结果处理标签、方向、比例尺/色键；当前仍是图位稿 |
| A14 | P2 / 全文 | 版本历史、待办、旧拟合误差混入正文；题名承诺超出有效结果 | 单独支持记录；重写题名、摘要、引言及结论，重排引用 | 事实审核后最终英文写作及具体期刊格式 |

## 可直接定位的原始与派生依据

- [实验室逐试样记录](/home/abcd/repos/work_hardening_paper/data/tensile_data/gr4b23271_lab_tensile_records.csv)：尺寸、标距、A/Z。
- [逐曲线提取汇总](/home/abcd/repos/work_hardening_paper/data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/gr4b23271_cold_deformation_tensile_summary.csv)：elastic_fit_range_MPa、elastic_modulus_GPa、elastic_intercept_MPa、elastic_fit_R2、elastic_fit_points、strain_at_Rp0.2、yield_status。
- [曲线口径的统一转写表](/home/abcd/repos/work_hardening_paper/tables/gr4b23271_mechanical_properties_curve_basis.csv)与[装配脚本](/home/abcd/repos/work_hardening_paper/tools/assemble_manuscript_mechanical_summary.py:68)：后者只装配已有结果，未重新拟合屈服。
- [旧图的截取程序](/home/abcd/repos/work_hardening_paper/tools/plot_manuscript_figures.py:315)及[负导数截零](/home/abcd/repos/work_hardening_paper/tools/plot_manuscript_figures.py:531)：用于识别旧图限制，不在本轮修改或运行。
- [GND 的 grainId 与 KAM 调用](/home/abcd/repos/work_hardening_paper/tools/mtex/generate_reference_style_fig14_gnd.m:96)、[Fig.6 分析脚本](/home/abcd/repos/work_hardening_paper/tools/mtex/generate_six_state_fig6_component_gallery.m)、[一阶锥面 Schmid 调用](/home/abcd/repos/work_hardening_paper/tools/mtex/generate_pyramidal_ca_schmid_distribution.m:15)。
- [原 Fig.6 统计](/home/abcd/repos/work_hardening_paper/results/mtex_fig6_component_gallery/component_summary.csv)、[原 GND 统计](/home/abcd/repos/work_hardening_paper/results/mtex_reference_style_fig14_gnd/reference_style_fig14_gnd_summary.csv)：方法限制改变解释范围，不代表本轮已重算数值。
- [材料命名与追溯](/home/abcd/repos/work_hardening_paper/docs/material_naming_and_traceability.md)：证书是钛锭资料，与本批棒材的对应尚待核对，不能直接引用其成分。

## KAM 实现与统计覆盖的直接核查

1. [Fig.6 生成器](/home/abcd/repos/work_hardening_paper/tools/mtex/generate_six_state_fig6_component_gallery.m:187) 先以 `boundary_detection_deg=2` 重构并赋 `ebsdFull.grainId`，再执行一阶、5° KAM。有限 Ti-Hex 像素取值在第 193–207 行。参数来自 [parameters.csv](/home/abcd/repos/work_hardening_paper/results/mtex_fig6_component_gallery/parameters.csv)。
2. [GND 生成器](/home/abcd/repos/work_hardening_paper/tools/mtex/generate_reference_style_fig14_gnd.m:94) 使用同样的 2°重构编号，再执行 KAM；第 105–106 行屏蔽非 Ti-Hex，相应有限数筛选位于第 121–123 行，`valid_pixel_count` 位于第 140 行。参数来自 [GND 参数表](/home/abcd/repos/work_hardening_paper/results/mtex_reference_style_fig14_gnd/reference_style_fig14_gnd_parameters.csv)。
3. 安装的 MTEX 6.1.1 [EBSDsquare/KAM.m](/home/abcd/MATLAB/toolboxes/mtex-6.1.1/EBSDAnalysis/@EBSDsquare/KAM.m:98) 第 98–104 行分别执行角度排除与跨 `grainId` 排除，两者同时生效；第 115 行除以有效邻点数。没有有效邻点时不生成有限 KAM。默认一阶方格邻域为四个正交近邻（同一文件第 41–44 行），不应笼统写成 8 邻点。

因此，原图 6c/9/14（现图 3c/6/11） 的 KAM/GND 不能解释为“不考虑 2°域边界、只用 5°截止计算的全场 KAM/GND”。它仍然能用于所采用条件下的域内梯度描述。不能仅由程序发现就断言当前所有 GND 结论失效，也不能认为现在补一句局限就已完成投稿级方法验证。

补充区别：用于 KAM 的这次 `calcGrains` 没有最少 5 像素筛选。最少 5 像素条件用于尺寸/Schmid 等相应保留区域统计，不应无依据推广为 KAM 只在不少于 5 像素的域中计算。必要时在方法中注明“KAM 所用域编号由未实施最少像素筛选的 2°重构获得”。

### 统计支持量（GND 模块，不挪作 原 Fig.6 模块的精确像素数）

分母为 [ODF 汇总](/home/abcd/repos/work_hardening_paper/results/mtex_odf_selected_phi2/odf_selected_phi2_summary.csv) 的 `valid_ti_hex_orientation_count`（同一 raw CTF 的已索引 Ti-Hex 取向数）；分子为 [GND 汇总](/home/abcd/repos/work_hardening_paper/results/mtex_reference_style_fig14_gnd/reference_style_fig14_gnd_summary.csv) 的 `valid_pixel_count`。全扫描为 600×600=360000 点。以下百分比只由这两个既有计数字段相除，不进行 EBSD 重处理。

| 名义直径/mm | 已索引 Ti-Hex 点数 | 有限 KAM-GND 点数 | 有限点/已索引 Ti-Hex | 有限点/全扫描 |
|---|---:|---:|---:|---:|
| 7.00 | 358203 | 357906 | 99.92% | 99.42% |
| 6.5 | 350441 | 320360 | 91.42% | 88.99% |
| 6.02 | 355266 | 346682 | 97.58% | 96.30% |
| 5.60 | 335389 | 285647 | 85.17% | 79.35% |
| 5.25 | 317224 | 236966 | 74.70% | 65.82% |
| 5.00 | 328099 | 257789 | 78.57% | 71.61% |

此处“未纳入有限 KAM”不等于“未索引”，后者已经在分母之前区分；也不能将所有失效像素归因为同一种物理缺陷。高变形量状态的有效支持范围不同，应纳入最终不确定性审查。

## CRSS 参考参数与滑移族对应

Xiong 作者全文第 2.3 节同时讨论一阶和二阶锥面，第 3.3 节将锥面 CRSS 固定为基面的三倍；表 4 的 756 MPa 仅以锥面族列示，没有分别标定两个阶次。本研究的一阶锥面来自 MTEX 的 `pyramidalCA` 定义，因此把 756 MPa 用于该族比较是本研究的参数对应假设，不能写为文献给出了一阶锥面的独立实测 CRSS。修订稿 4.2 已明确这一点。[Xiong 作者全文](https://arxiv.org/pdf/2003.01682)

## 与最接近论文的关系

Meng 等已经将无热处理纯钛旋锻的强度变化与晶界、位错和织构联系起来。因此本稿的研究价值应落在实际初态、变形范围和组织演化差异上，并用重复取样及局部证据检验，不能以增加表征种类替代知识增量。[原始论文](https://www.sciencedirect.com/science/article/pii/S0925838820345850)

## 配套记录的纠正

本轮同步纠正前次 MD 综合评估中“旋锻态约 87%–90%”的过宽概括。正确值为：五个旋锻状态的 LAGB 长度分数依次为 76.16%、53.48%、85.16%、90.19%、87.18%；87.18%–90.19% 只对应最高两档减缩率，不代表全部旋锻状态。正文表内数值原本正确。

参考文献目录中关于本项目 O/Fe 含量的旧措辞也需与证书追溯保持一致；本轮将其标为未确认对应的旧记录，避免被当作本批实测成分再次引用。

## 核对结论的范围

本轮只读核对显示，六状态三个曲线均值的 18 项聚合结果与既有提取表一致；已有组织分布字段核对通过。它们验证的是汇总与转写，不验证原始仪器通道、屈服算法、EBSD 重构的全部适用性，也不能弥补独立实验重复不足。

最终文档一致性检查 25/25 项通过，包括：原稿 SHA-256 未变化；当前四张表的 24 行状态数据与原稿一致；旧模型六行数值完整保留；18 个原图位在正文与支持记录中均可追溯，其中 17 个既有图片文件存在，TEM 为唯一未生成的计划图；正文图表及七条参考文献编号连续，当前交付文档的本地链接可解析；六状态有效像素比例可由既有计数字段复现。正文使用 12 张既有图件和 1 个 TEM 图位，另 5 张旧图转入支持记录。本轮未制作新图，也未修改原始数据。

检查脚本及完整结果分别保存于[本轮文档检查脚本](/home/abcd/repos/work_hardening_paper/.codex_tmp/2026-09-08-paper-review/verify_revision.py)和[检查结果](/home/abcd/repos/work_hardening_paper/.codex_tmp/2026-09-08-paper-review/revision_verification.json)。这些是文档一致性检查，不是对实验或机制结论的独立验证。

人工处理顺序与各模块完成标志见[模块审核清单](/home/abcd/repos/work_hardening_paper/docs/2026-09-08-manuscript-module-review.md)。
