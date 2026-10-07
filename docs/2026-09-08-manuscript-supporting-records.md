# 历史图件、强化模型与编写记录

日期：2026-09-08。本文档用于项目追溯，不属于修订稿正文，也不自动作为投稿补充材料。当前正文见[修订稿](/home/abcd/repos/work_hardening_paper/manuscript_revision_2026-09-08.md)，最新方法审查见[证据审核](/home/abcd/repos/work_hardening_paper/docs/2026-09-08-manuscript-evidence-audit.md)。以下旧文中的图号、参考编号、参数来源及待办状态均按原稿保存；特别是KAM方法应以新稿披露的2°域内限制为准。

## A. 转入记录的旧图位

旧图2的强度采用实验室记录；旧图3/4的截止与求导处理不适于当前定量硬化结论；旧图15/16使用旧强度基准及未独立确认的模型参数。文件完整保留，未重绘、改标签或生成替代图。

> **图片占位｜图 2：六状态力学性能汇总图（旧版图位）**  
> 路径：[ta4_mechanical_property_summary.svg](/home/abcd/repos/work_hardening_paper/figures/ta4_mechanical_property_summary.svg)  
> 图注：既有实验室记录口径的强度与断后延性汇总。当前正文 Rp0.2、Rm 采用表 2a 的曲线口径，A、Z 采用表 2b；旧图中的实验室强度数值不作为表 2a 的图示结果。  
> 来源：项目已有论文图；本次保留占位，不重绘。当前强度的定量依据为表 2a。


> **图片占位｜图 3：六状态代表性真实应力–应变曲线**  
> 路径：[ta4_true_stress_strain_representative_curves.svg](/home/abcd/repos/work_hardening_paper/figures/ta4_true_stress_strain_representative_curves.svg)  
> 图注：工程量换算的真实应力–应变图；旧图采用最大换算真应力对应应变作为截取上限，超过工程应力峰值的部分不解释为局部真实应力–应变。  
> 来源：项目已有论文图，保留原文件；6.48 mm 状态按本文名义分组对应 Y-6.5。


> **图片占位｜图 4：六状态代表性加工硬化率曲线**  
> 路径：[ta4_work_hardening_rate_representative_curves.svg](/home/abcd/repos/work_hardening_paper/figures/ta4_work_hardening_rate_representative_curves.svg)  
> 图注：经重采样、平滑与中心差分得到的加工硬化率显示值，横轴为总真实应变；负导数被截为零，曲线终点未统一对应工程应力峰值。  
> 来源：项目已有论文图，本次仅保留图位，不据其给出硬化阶段或回复参数。


> **图片占位｜图 15：文献参数版的基体、固溶、HAGB 与位错强度分解**  
> 路径：[reference_three_factor_strengthening_montage.png](/home/abcd/repos/work_hardening_paper/results/reference_three_factor_strengthening/figures/reference_three_factor_strengthening_montage.png)  
> 图注：基体、固溶、HAGB 和 KAM-GND 位错项的文献参数估算及实验屈服强度比较；各状态采用相同文献参数，未使用六状态屈服强度进行参数拟合。  
> 来源：项目已有强化贡献图，作为图 16 校准结果的参照。


> **图片占位｜图 16：共用参数校准后的基体、固溶、HAGB 与位错强度分解**  
> 路径：[reference_three_factor_strengthening_calibrated_montage.png](/home/abcd/repos/work_hardening_paper/results/reference_three_factor_strengthening/figures/calibrated/reference_three_factor_strengthening_calibrated_montage.png)  
> 图注：两个共用参数针对六状态屈服强度校准后的贡献分解，基体项为 76.75 MPa，有效密度尺度比为 3.437，固溶及 HAGB 参数保持不变。图示为描述性拟合结果。  
> 来源：项目图件 PPT 第 13 页。


## B. 修订前强化模型原文（历史资料，不作为当前结果）

### 4.4 强化贡献与强塑性关系

上述分析表明，晶粒尺度、取向和位错结构分别通过晶界约束、分切条件及滑移阻力影响力学响应。为进一步判断旋锻引起的强度增加主要来自何种组织变化，需要将各项的绝对强度与相对于未变形态的强化增量区分开来。

> **本节数据注记（不属于论文正文）：** 图 15、16 和表 4 保留项目既有模型，比较基准为旧实验室屈服强度，与表 2a 的曲线复算口径不同。本节数值用于搭建贡献讨论，所列误差不用于验证当前曲线结果。模型包含基体、固溶、HAGB 和位错项，固定 Taylor 因子，尚未独立量化取向贡献。

在同批材料组成保持不变的假设下，基体及固溶项作为共同强度基准，HAGB 晶界项与位错项描述所采用组织参数对强度的影响。采用与文献 [5] 相同的叠加形式：

\[
\sigma_y=\sigma_0+\sigma_{\mathrm{SS}}
+\sigma_{\mathrm{HP}}+\sigma_{\mathrm{dis}},\tag{3}
\]

\[
\sigma_{\mathrm{HP}}=\frac{k}{\sqrt d},\qquad
\sigma_{\mathrm{dis}}=\alpha M G b\sqrt{\rho_{\mathrm{GND}}}.\tag{4}
\]

式中，σ₀ 为基体应力项，σSS 为固溶项，σHP 为 HAGB 晶界项，σdis 为基于 KAM-GND 的位错项。文献参数版取 σ₀=172.5 MPa、σSS=258.11 MPa、k=330 MPa·μm¹ᐟ²、α=0.2、M=5.0、G=45.6 GPa、b=0.295 nm；d 使用 15°HAGB 晶粒 ECD 的数量均值，单位为 μm。σ₀、α、M 及 b 参考 Xia 等采用的参数[5]，k 和 G 参考纯 α-Ti 的晶界强化分析[7]，固溶项沿用基于间隙溶质模型的已有估算[6]。M 在各状态中固定为 5.0，未根据本研究取向数据逐状态计算。

图 15 和表 4 给出文献参数下的强化分解。HAGB 项为 132.0–138.9 MPa，与晶粒数量平均 ECD 的有限变化相对应；位错项则由未变形态的 76.9 MPa 增至最终态的 246.6 MPa。固溶项在各状态中保持 258.11 MPa，不提供随减缩率变化的增量。因此，在这一固定参数模型中，位错项承担了主要的状态差异。

> **图片占位｜图 15：文献参数版的基体、固溶、HAGB 与位错强度分解**  
> 路径：[reference_three_factor_strengthening_montage.png](/home/abcd/repos/work_hardening_paper/results/reference_three_factor_strengthening/figures/reference_three_factor_strengthening_montage.png)  
> 图注：基体、固溶、HAGB 和 KAM-GND 位错项的文献参数估算及实验屈服强度比较；各状态采用相同文献参数，未使用六状态屈服强度进行参数拟合。  
> 来源：项目已有强化贡献图，作为图 16 校准结果的参照。

**表 4 已有强化模型与旧实验室屈服强度基准的比较（单位：MPa）**

| 直径 / mm | 旧实验室 Rp0.2 | HAGB 项 | 文献参数位错项 | 文献参数总强度 | 校准位错项 | 校准总强度 |
|---|---:|---:|---:|---:|---:|---:|
| 7.00 | 580.5 | 138.9 | 76.9 | 646.4 | 142.5 | 616.2 |
| 6.5 | 845.0 | 132.0 | 232.0 | 794.7 | 430.2 | 897.1 |
| 6.02 | 900.0 | 133.2 | 203.0 | 766.8 | 376.4 | 844.4 |
| 5.60 | 927.0 | 132.9 | 242.9 | 806.4 | 450.4 | 918.1 |
| 5.25 | 982.5 | 138.7 | 248.2 | 817.5 | 460.1 | 933.7 |
| 5.00 | 980.0 | 138.0 | 246.6 | 815.2 | 457.2 | 930.1 |

表注：文献参数版基体项为 172.5 MPa，校准版为 76.75 MPa；两版固溶项均为 258.11 MPa。表中为既有计算结果的舍入值，各项之和可能存在 0.1 MPa 量级的显示舍入差。

文献参数模型高估未变形态、总体低估旋锻态相对于旧实验室记录的强度，平均绝对百分比误差为 13.12%，最大绝对相对偏差为 16.81%。组合文献参数的适用性及 KAM-GND 对位错结构的有限表征均可能影响绝对强度估算。尤其是固溶项所用氧当量 0.3484 wt.% 超出来源模型的 0.14–0.32 wt.% 验证区间，且未单独处理 Fe 的贡献[6]；该成分输入与本批棒材的对应关系尚未确认。因此，固溶项只能作为既有估算保留，不能视为本批材料已验证的强化值。

图 16 给出采用两个共用参数校准后的结果。保持固溶项与 HAGB 参数不变，对六状态采用共同参数校准，得到有效基体应力 76.75 MPa 和位错强化倍率 1.854，对应有效密度尺度比 3.437。相对于同一旧实验室强度基准，平均绝对百分比误差降至 4.92%，最大绝对相对偏差为 6.17%；位错项由 142.5 MPa 增至 457.2 MPa，仍承担主要的模型状态差异。校准仅提高了对参与拟合数据的描述能力，有效密度尺度比不是实测总位错密度与 GND 的比值，有效基体应力也不直接等同于晶格摩擦应力。

> **图片占位｜图 16：共用参数校准后的基体、固溶、HAGB 与位错强度分解**  
> 路径：[reference_three_factor_strengthening_calibrated_montage.png](/home/abcd/repos/work_hardening_paper/results/reference_three_factor_strengthening/figures/calibrated/reference_three_factor_strengthening_calibrated_montage.png)  
> 图注：两个共用参数针对六状态屈服强度校准后的贡献分解，基体项为 76.75 MPa，有效密度尺度比为 3.437，固溶及 HAGB 参数保持不变。图示为描述性拟合结果。  
> 来源：项目图件 PPT 第 13 页。

在上述共同固定基体项和固溶项的模型中，相对于未变形态的强化增量为：

\[
\Delta\sigma_y=\Delta\sigma_{\mathrm{HP}}
+\Delta\sigma_{\mathrm{dis}}.\tag{5}
\]

两套参数下 HAGB 项变化均较小，位错相关项均构成主要模型增量。Xia 等的分解则显示，绝对位错项较大与所比较冷轧状态间的晶界增量较大可以同时成立[5]。两者说明，“贡献最大”必须明确比较的是绝对强度还是变形引起的增量。本研究中，文献参数版的恒定固溶项仍高于各状态位错项，因而位错主导的判断应限定于模型中的旋锻强化增量。

从组织与性能的对应关系看，HAGB 晶粒尺寸变化有限，取向演化调整了滑移几何条件，而晶内低角度亚结构和位错相关取向梯度整体增加。结合强化分解，现有结果支持位错相关结构是旋锻强度提高的主要组织来源之一，其中位错项在当前模型中承担最大增量。由于模型采用固定 M，取向影响未被独立分离，校准参数也可能吸收未计位错和文献参数差异，尚不能据此给出晶粒、取向和位错三者的独立贡献百分比或完整定量排序。

强度增加同时伴随断后伸长率由 26.75% 降至 8.00%、断面收缩率由 44.50% 降至 31.00%（表 2b）。位错相互作用和晶内边界约束提高滑移阻力，也可能限制后续变形协调，从而与这一强塑性变化相联系[1,5]。但较高的预变形强度并不直接意味着后续拉伸中的加工硬化率更高，屈服强度分解也不预测断后延性。当前结果支持旋锻强化与后续延性降低的总体对应关系；其加工硬化阶段及储存–回复平衡仍需可靠的拉伸硬化区间和相应组织证据加以判定。

## C. 旧编号参考文献（仅供本文件旧文对应）

[1] Wang M, Wang Y, Huang A, Gao L, Li Y, Huang C. Promising tensile and fatigue properties of commercially pure titanium processed by rotary swaging and annealing treatment. *Materials*, 2018, 11(11): 2261. [DOI: 10.3390/ma11112261](https://doi.org/10.3390/ma11112261).

[2] Meng A, Chen X, Nie J, Gu L, Mao Q, Zhao Y. Microstructure evolution and mechanical properties of commercial pure titanium subjected to rotary swaging. *Journal of Alloys and Compounds*, 2021, 859: 158222. [DOI: 10.1016/j.jallcom.2020.158222](https://doi.org/10.1016/j.jallcom.2020.158222).

[3] Xiong Y, Karamched P, Nguyen C-T, Collins DM, Magazzeni CM, Tarleton E, Wilkinson AJ. Cold creep of titanium: Analysis of stress relaxation using synchrotron diffraction and crystal plasticity simulations. *Acta Materialia*, 2020, 199: 561–577. [DOI: 10.1016/j.actamat.2020.08.010](https://doi.org/10.1016/j.actamat.2020.08.010).

[4] El-Dasher BS, Adams BL, Rollett AD. Viewpoint: experimental recovery of geometrically necessary dislocation density in polycrystals. *Scripta Materialia*, 2003, 48(2): 141–145. [DOI: 10.1016/S1359-6462(02)00340-8](https://doi.org/10.1016/S1359-6462(02)00340-8).

[5] Xia Y, Yu G, Yoon B-H, Cao P, Yi Y, Yu Z. Numerical simulation of cold working and microstructure–property regulation of TA16 titanium alloy. *Materials Science and Engineering: A*, 2026, 975: 150869. [DOI: 10.1016/j.msea.2026.150869](https://doi.org/10.1016/j.msea.2026.150869).

[6] Won JW, Park CH, Lee T, Lee CS. Integrated constitutive model for flow behavior of pure Titanium considering interstitial solute concentration. *Metals and Materials International*, 2014, 20(6): 1017–1025. [DOI: 10.1007/s12540-014-6004-8](https://doi.org/10.1007/s12540-014-6004-8).

[7] Luo P. Analysis of microstructure and its effect on yield strength of pure alpha-titanium consolidated by equal channel angular pressing. *Materials Transactions*, 2018, 59(7): 1161–1165. [DOI: 10.2320/matertrans.M2018033](https://doi.org/10.2320/matertrans.M2018033).

## D. 修订前编写注记原文（历史资料）

## 编写注记（不属于论文正文）

### A. 本稿图件安排

本稿按出现顺序预留图 1–16。其中 17 个已有图片文件对应图 1–12、14–16，图 6a、6b、6c 分别保留 IPF、GB、KAM 原图；图 13 仅为 TEM 名称与计划路径，没有对应文件。图 11、12 为柱面〈a〉和锥面〈c+a〉施密特因子分布，图 14 为 GND，图 15、16 为已有强化模型图。本轮没有重绘、拼接或生成图片。

PPT 来源：[gr4b23271_reference_figures_fig6_fig15.pptx](/home/abcd/repos/work_hardening_paper/results/presentations/gr4b23271_reference_figures_fig6_fig15.pptx)。其中第 8 页是 PF/IPF，第 11 页是初末态四滑移族空间图；本稿按指定结构选用已有 ODF 截面及两张 Schmid 分布图，未强行加入上述两页。

### B. 资料对应与已知待核对项

1. **实验参数的来源层级。** 厂家、轴向室温拉伸、INSTRON 5582、GB/T 228.1-2010、室温旋锻及无中间退火、Kroll 腐蚀、机械与电解抛光等表述沿用 [原稿实验部分](/home/abcd/repos/work_hardening_paper/manuscript_draft.md:25)。本轮以实验记录另行确认了 L₀=Lₑ=50 mm、逐试样实测直径和重复数；以六个 raw CTF 头信息确认了 20 kV、70°、600×600 和 0.5 μm。拉伸速率、旋锻道次与速度、温升与润滑、SEM/EBSD 型号及电解抛光参数仍无明确记录，未补写推测值。
2. **材料名称及资料追溯。** 本稿采用“Grade 4 工业纯钛”，英文为“Grade 4 commercially pure titanium”，需要缩写时使用“Grade 4 CP-Ti”。内部生产编号仅用于项目资料追溯，不进入论文标题、摘要、正文、图注或试样名称；M-7.00 与 Y-d 仅表示处理状态和名义直径。旧稿及文件名中的材料标记不直接用于正式命名，原始数据和路径保持原样。当前棒材与材料证明的对应关系仍待核对，证书成分不直接作为本批试样实测成分，也不据现有名称宣称已确认特定材料标准的符合性。内部标识、历史差异、证书信息及来源记录见 [材料命名与内部追溯记录](/home/abcd/repos/work_hardening_paper/docs/material_naming_and_traceability.md)。
3. **名义分组与旧图标记。** 论文名称统一为 6.5 mm，以名义直径计算 ψ=13.78%；曲线、EBSD 和部分旧图中的 6.48 mm／14.31% 作为同一状态的历史标记保留。拉伸记录中该组两根实测直径为 6.52、6.51 mm。未改动原始文件名、CTF、应力应变字段或既有图内文字。金相图 5 的方向箭头仍标为 RD，与正文 AD 的映射尚需核对，本稿仅描述其定向拉长。
4. **曲线筛选与断后记录。** 当前摘要、表 2a 及力学正文的 Rp0.2、Rm 均采用人工排除后的曲线汇总。6.02 Y-2 的原始屈服状态为 `not_reliable_no_elastic_segment_or_offset_crossing`，按已有人工排除决定剔除整条曲线的强度和峰值应变；该组有效曲线 n=1，标准差记为“—”。A、Z 采用实验室断后记录，六状态各 n=2。表 2b 中 Y-6.02 的 A=10.00±0.00 表示两条同精度记录相同，不是单值汇总的零标准差。新的派生表保留所有来源、样本名称和排除原因，未覆盖旧表。
5. **旧力学图的适用范围。** 图 2 仍为实验室强度口径，不能充当表 2a 新曲线口径的图示；本次按“不生成新图”的约定保留占位。图 3、4 对应的 [绘图脚本](/home/abcd/repos/work_hardening_paper/tools/plot_manuscript_figures.py:306) 以最大换算真应力对应应变截取，并在硬化率中对负导数截零，因此不能用整幅曲线确定颈缩前区间或硬化阶段。本轮定量表述采用 εe,Rm、强度及断后延性，没有重算或修改这些图。
6. **EBSD 统计口径。** 图 6–9 使用 Fig.6 图组的 2°检测下限，不混入另一套 1°下限的 LAGB 分数。HAGB 尺寸分布和既有强化模型采用数量口径，表 3 另列面积加权中位数。Schmid 图文件的 `grain` 字段实际对应 2°重构取向域。GND 使用其自身模块的 KAM 数据，未将 Fig.6 模块未舍入的 KAM 数值代入替换。
7. **强化模型与本轮力学口径的关系。** 第 4.4 节、表 4 和图 15、16 保留之前的模型及旧实验室屈服强度基准。本轮没有以表 2a 的新曲线强度重算模型或重新拟合，因此旧 MAPE=4.92%、最大偏差 6.17% 不用于宣称对当前曲线口径的拟合一致性。文献 CRSS、固溶外推和有效密度比 3.437 的解释边界继续保留，后续讨论修订需与选定力学口径对应。

8. **三类组织因素与既有模型项的区别。** 本轮讨论对象为晶粒/晶界、取向、位错；既有“three_factor”图中的三个项实际为固溶、HAGB、位错，另有基体项。固定 M=5 的图不能作为独立取向贡献结果。最终比较针对相对于 M-7.00 的强化增量，未给出三项贡献百分比，也未把当前模型描述性结果写成已验证的排序。
9. **TEM 的位置与状态。** 方法预留并入 2.3 节，结果位于 ODF/Schmid 后的 3.4 节，后接 3.5 节 GND。图 13 是唯一尚无文件的计划图位；取样状态、数量及条件均未锁定。TEM 重点补充局部晶体学关系与位错组态，宏观取向演化仍以 EBSD 统计为依据。若“转变”具体指滑移族主导关系变化，还需位错类型、滑移面及代表性比较；不能将预期观察预写为结果。方法参考为文献 [5] 的 3.2 节/Fig.8 及 4.2 节/Fig.12，其 TA16 观察结论不移用于本材料。

### C. 关键数值来源索引

| 本稿内容 | 已有来源文件 |
|---|---|
| 本轮统一分组、曲线强度与断后延性派生表 | [gr4b23271_mechanical_properties_curve_basis.csv](/home/abcd/repos/work_hardening_paper/tables/gr4b23271_mechanical_properties_curve_basis.csv) |
| 本轮派生表的复现脚本 | [assemble_manuscript_mechanical_summary.py](/home/abcd/repos/work_hardening_paper/tools/assemble_manuscript_mechanical_summary.py) |
| 拉伸逐试样尺寸、标距与断后记录 | [gr4b23271_lab_tensile_records.csv](/home/abcd/repos/work_hardening_paper/data/tensile_data/gr4b23271_lab_tensile_records.csv) |
| 实验室转写工作簿及转写说明 | [lab_records.xlsx](/home/abcd/repos/work_hardening_paper/data/tensile_data/lab_records.xlsx) |
| 材料质量证明（与棒材批次对应尚待核对） | [gr4b23271_grade4_quality_certificate.png](/home/abcd/repos/work_hardening_paper/data/material_certificate/gr4b23271_grade4_quality_certificate.png) |
| EBSD 采集条件（六状态均已检查，此列以初态为例） | [ebsd_sample_7_map_15.ctf](/home/abcd/repos/work_hardening_paper/data/ebsd_kpl_250221_7_df/scans/d7/ebsd_sample_7_map_15.ctf:5) |
| 断后 A/Z 均值与旧模型所用实验室强度基准 | [gr4b23271_lab_tensile_by_diameter.csv](/home/abcd/repos/work_hardening_paper/data/tensile_data/gr4b23271_lab_tensile_by_diameter.csv) |
| 曲线复算、人工排除后的状态结果 | [gr4b23271_tensile_final_by_diameter_user_exclude.csv](/home/abcd/repos/work_hardening_paper/data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/gr4b23271_tensile_final_by_diameter_user_exclude.csv) |
| 代表性曲线选择与原始文件对应 | [ta4_representative_tensile_curves.csv](/home/abcd/repos/work_hardening_paper/figures/ta4_representative_tensile_curves.csv) |
| 15°HAGB 晶粒尺寸 | [hagb_grain_size_summary_used.csv](/home/abcd/repos/work_hardening_paper/results/mtex_fig6_component_gallery/07_hagb_grain_size/hagb_grain_size_summary_used.csv) |
| 2°取向域、LAGB 分数、KAM | [component_summary.csv](/home/abcd/repos/work_hardening_paper/results/mtex_fig6_component_gallery/component_summary.csv) |
| ODF 全空间最大值及参数 | [odf_selected_phi2_summary.csv](/home/abcd/repos/work_hardening_paper/results/mtex_odf_selected_phi2/odf_selected_phi2_summary.csv) |
| 柱面〈a〉分布 | [prismatic_a_schmid_summary.csv](/home/abcd/repos/work_hardening_paper/results/mtex_prismatic_schmid_distribution/prismatic_a_schmid_summary.csv) |
| 锥面〈c+a〉分布 | [pyramidal_ca_schmid_summary.csv](/home/abcd/repos/work_hardening_paper/results/mtex_pyramidal_ca_schmid_distribution/pyramidal_ca_schmid_summary.csv) |
| 固定 CRSS 下的逐取向域竞争 | [ca_vs_a_preference_summary.csv](/home/abcd/repos/work_hardening_paper/results/mtex_ca_vs_a_preference/ca_vs_a_preference_summary.csv) |
| GND 均值及估算参数 | [reference_style_fig14_gnd_summary.csv](/home/abcd/repos/work_hardening_paper/results/mtex_reference_style_fig14_gnd/reference_style_fig14_gnd_summary.csv) |
| 文献参数与校准强化贡献 | [reference_three_factor_strengthening.csv](/home/abcd/repos/work_hardening_paper/results/reference_three_factor_strengthening/reference_three_factor_strengthening.csv) |
| 强化模型参数及误差 | [reference_three_factor_strengthening.json](/home/abcd/repos/work_hardening_paper/results/reference_three_factor_strengthening/reference_three_factor_strengthening.json) |

### D. 各模块的写作任务

| 模块 | 需要回答的问题 | 采用材料或预留内容 |
|---|---|---|
| 研究背景 | 为什么研究旋锻冷变形，现有研究还不能解释什么 | 强塑性调控需求及晶粒、取向、位错相对作用的研究问题 |
| 实验内容 | 对什么材料做了什么，比较条件是否一致 | 六状态分组、拉伸、金相、EBSD；TEM 明确列为计划方法 |
| 3.1 性能 | 旋锻后强度、延性与后续拉伸响应如何变化 | 曲线、强度和断后延性共同组成拉伸性能；后续拉伸响应单列，保留旧硬化图的适用范围 |
| 3.2 晶粒与亚结构 | 是否出现晶粒尺度变化，晶内边界如何演化 | 金相和 IPF/GB/KAM 展示形貌，再以尺寸、晶界角度及 KAM 分布比较；相关图件合为连续论述 |
| 3.3 取向 | 晶体取向及轴向滑移几何条件如何变化 | ODF 与滑移几何条件两个单元；柱面与锥面 Schmid 分布并列比较 |
| 3.4 TEM 预留 | 局部取向关系和位错结构是否支持相关解释 | 衍射标定、局部亚结构及必要的位错类型分析，全部待补 |
| 3.5 GND | 可测取向梯度对应的位错密度尺度如何变化 | 已有 KAM-GND 图；与 KAM 非独立、与总位错密度有区别 |
| 4.1 晶粒与晶界 | HAGB 晶粒未连续细化时，晶界结构如何参与强化 | 区分高角度晶粒尺度与低角度亚结构，并与文献组织演化比较 |
| 4.2 取向与滑移 | 高施密特因子为何能与较高强度共存 | 先提出结果中的对应关系，再以分切条件和 CRSS 解释；TEM 关联内容待补 |
| 4.3 位错与阻力 | 位错相关结构如何提高阻力，局部非单调变化如何理解 | 引入 KAM-GND 与位错相互作用，说明取样、密度估算及后续硬化解释的范围 |
| 4.4 强化与强塑性 | 模型中哪一项主导增量，这与延性变化有何联系 | 集中比较两套既有分解；区分绝对强度和增量，保留旧强度基准与取向未独立量化的说明 |
| 结论 | 哪些结果已经成立，主要强化来源能否确定 | 保留已测结论；位错贡献最大的最终措辞附验证条件 |

### E. 参考论文的章节体例与本稿对应

本轮对照 [Xia 等 TA16 冷加工论文](/home/abcd/repos/work_hardening_paper/references/xia_2026_ta16_cold_rolling_microstructure_property.pdf) 的第 3、4 节，采用其按组织或性能对象归并结果、围绕变形和强化机制展开讨论的组织方式。其 3.2 节将 EBSD 形貌、尺寸、晶界及取向统计与 TEM 观察连续安排，3.3 节集中报告力学性能；4.2 节再结合 Schmid/CRSS 和 TEM 讨论滑移，4.3 节由组织变化引入强化关系并比较各项贡献。本稿沿用用户确定的性能在先、组织在后的顺序。

| 内容 | 本稿位置 | 写作层级 |
|---|---|---|
| 材料来源、尺寸、旋锻及热处理状态 | 2.1 | 交代材料和加工条件 |
| 拉伸仪器、取样方向、标距、重复数及基本筛选 | 2.2 | 交代试验方法及有效样本 |
| 金相、EBSD 制样与采集条件，TEM 待做说明 | 2.3 | 交代显微组织表征方法 |
| 拉伸曲线、强度、断后延性及后续响应 | 3.1 | 按拉伸性能与后续响应组织，图表服务于状态比较 |
| 金相、IPF/GB/KAM 和组织统计 | 3.2 | 由形貌观察到定量比较，取消每张图对应一个小节的安排 |
| ODF 与两类施密特因子分布 | 3.3 | 先比较取向密度，再并列比较两类滑移族的几何条件 |
| TEM 预留及 GND 分布 | 3.4、3.5 | 保留计划观察对象，报告已有密度统计；机制推论归入讨论 |
| 晶粒尺度与晶内低角度细分的区别 | 4.1 | 从 HAGB 尺寸与强度不同步的结果出发，讨论边界作用 |
| Schmid 定义、载荷假设及 CRSS 比较 | 4.2 | 先说明高 SF 与高强度共存的问题，再引入公式和文献解释 |
| KAM-GND 关系及位错结构作用 | 4.3 | 从亚结构增加引入密度估算，再讨论滑移阻力与局部差异 |
| 强化分解、参数及强塑性关系 | 4.4 | 模型估算、参数校准和物理归纳连续展开，避免把校准过程写成独立研究主题 |
| 名义减缩率换算、曲线处理及详细统计规则 | 编写注记 F | 保存计算记录，不占据实验正文 |

Results 中保留图示现象、重点数值比较及观察层面的归纳；滑移阻力、晶界作用、位错强化以及组织与强度局部不一致的解释集中于 Discussion。各讨论小节由本研究结果提出具体问题，再引入相关物理关系和文献，最后限定可得结论。Schmid、GND 及强化公式均用于解释相应问题，不作为小节的直接开场。

本文未开展参考论文中的有限元模拟，不设置相应方法小节；TEM 继续位于取向与 GND 之间，实际观察及机制讨论均待补。参考论文的连续晶粒细化、孪生和滑移转变结论不直接移用于本研究。其强化分解还区分了绝对位错项较大与所比较状态间晶界增量较大，本稿据此明确“位错贡献最大”所指的比较尺度，避免将参考论文结论直接作为本材料结论。

### F. 数据处理与统计记录（不属于论文正文）

#### F.1 名义减缩率与尺寸口径

以 D₀=7.00 mm 为初始名义直径，D₁ 为各状态名义分组直径，减缩率按下式计算：

\[
\psi=\left[1-\left(\frac{D_1}{D_0}\right)^2\right]\times100\%.\tag{S1}
\]

名义尺寸只用于分组和变形量比较，应力仍采用各根试样的实测原始截面积。6.48 mm／14.31% 为旧文件与图件的历史标签，论文分组统一为 6.5 mm／13.78%，未对原曲线应力、应变或试样尺寸重新缩放。

#### F.2 基础应力–应变换算及既有曲线处理

均匀变形阶段，工程应力 σe 和工程应变 εe 按下式换算为真实应力 σt 与真实应变 εt：

\[
\sigma_t=\sigma_e(1+\varepsilon_e),\qquad
\varepsilon_t=\ln(1+\varepsilon_e).\tag{S2}
\]

加工硬化率定义为 θ=dσt/dεt。已有曲线采用总真实应变，重采样间隔为 0.00025，采用约 0.0025 应变宽度的移动平均及约 0.002 应变半宽的中心差分。工程应力达到最大值之前的均匀变形区间是式（S2）及加工硬化讨论的适用范围，颈缩后按相同公式继续换算的量不解释为局部真实应力与应变。

既有派生曲线以最大换算真应力对应应变为截取上限，且加工硬化率的负导数被截为零，其完整显示范围并不严格对应颈缩前区间。因此，本研究采用工程应力峰值对应的记录应变 εe,Rm、强度水平及断后延性描述后续拉伸响应，不据旧导数图识别硬化阶段转折、回复参数或饱和应力。εe,Rm 未扣除弹性应变，且未完成跨试样应变零点的统一复核，故不将其标记为标准塑性延伸率 Ag。

#### F.3 曲线筛选及断后数据口径

6.02 mm 状态的第 2 条曲线在接近零工程应变时已记录约 815.5 MPa 的应力，缺少能够可靠识别的弹性段或 0.2% 偏移交点，因此将该曲线从强度与曲线应变统计中排除。筛选后共纳入 11 条曲线：Y-6.02 仅保留第 1 条曲线，其余五种状态各保留两条。Rp0.2、Rm 及 Rm 对应工程应变均采用同一套纳入规则，避免不同指标使用不同曲线集合。

对于具有两条有效曲线的状态，报告算术平均值与样本标准差；Y-6.02 的曲线指标报告单条有效结果，标准差记为“—”，不将单值汇总中的零值解释为无离散性。强度结果统一采用曲线复算口径，不与实验室记录中的屈服强度混合求平均。名义 Y-6.5 状态在源曲线中的名称仍为“6.48 Y-1”和“6.48 Y-2”，仅调整论文分组标签。

断后伸长率 A 和断面收缩率 Z 采用实验室断后测量记录，分别对应 A=(Lᵤ−L₀)/L₀×100% 和 Z=(S₀−Sᵤ)/S₀×100%，其中 Lᵤ 为断后标距，S₀ 为原始截面积，Sᵤ 为断口处最小截面积。曲线起始段的记录问题不作为剔除断后测量值的依据，因此各状态的 A、Z 均保留两条实验室记录。曲线末点应变和 Rm 对应工程应变不替代 A；亦不根据曲线应变反算 Z。

各状态的代表性曲线为既有选择表对应的第 1 条曲线；组均值与其对应峰值分别统计，不能据代表曲线峰值要求组均值逐状态具有相同变化。

#### F.4 EBSD 与取向分布的详细统计设置

晶粒与 2°取向域的等效圆直径按 ECD=2√(S/π) 计算，其中 S 为重构区域面积；尺寸直方图采用数量频率，并另列面积加权中位数。原始与去噪数据不作为独立重复，晶粒或像素数量也不等同于独立试样数。

KAM 空间图采用统一 0°–5°色标，分布按有效 Ti-Hex 像素计数，分箱宽度为 0.2°。ODF 各状态采用统一 MRD 色标，全取向空间最大值由 1°分辨率寻峰获得，不与单截面峰值混用。Schmid 分箱宽度为 0.02，m=0.4 为统计参考线；原始数据作为主结果，去噪结果用于处理敏感性比较。GND 采用其对应模块的 KAM 与有效像素统计结果，不以另一模块的未舍入 KAM 数值替换。

### G. Results 篇幅与图文深度

本轮保留现有章节结构，增加各组图的空间与分布特征、重点状态比较及观察层面的解释。重点扩充图 6a–c、晶粒与晶界统计、ODF、两类 Schmid 和 GND；图 3、4 因既有处理限制维持有限解释，TEM 继续留空。篇幅不按每张图等额分配，也不将参考论文的英文词数换算成中文硬性字数。逐图审查、修订前后正文篇幅及新增数值的字段和权重记录见 [Results 篇幅与图文深度审查](/home/abcd/repos/work_hardening_paper/docs/2026-09-06-results-length-and-figure-depth-review.md)。

新增晶界区间比例来自既有长度频率分箱，KAM 与 GND 来自各自模块的像素统计；Schmid 高值区间采用取向域面积权重。晶界 4°–10°、Schmid 0.48–0.50 及 GND≥4×10¹⁴ m⁻² 均用于描述分布差异，不作为机制启动或饱和阈值。中位数与百分位数读取已有统计字段，未重新进行 EBSD 重构或概率分布拟合。

