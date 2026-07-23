# Gr4B23271 冷变形 EBSD 综合分析

本目录由统一的 MTEX 流程从六组原始 CTF 及其配对去噪数据生成。图中水平方向为棒材轴向 AD（试样 x 方向），竖直方向为面内横向/径向 TD/RD（试样 y 方向），扫描法向为 ND（试样 z 方向）。

## 数据作用与分析参数

原始 CTF 是定量结论的主要来源；去噪 CTF 仅用于逐状态敏感性比较，不作为独立重复。主要晶粒重构采用 2° 检测下限、15° 高角度晶界划分阈值及 minPixel=5，不进行晶界平滑。晶界统计按物理线段长度加权，同时保留 0.5°、1°、2° 和 5° 检测下限敏感性结果。

KAM 采用一阶邻域和 5° 排除阈值作为主要描述，同时输出一阶/二阶邻域与 2°/5° 阈值组合的敏感性结果。GROD 轴统计仅纳入 GROD 角不低于 5° 的像素，并以反极性二阶矩张量汇总轴向；KAM、GROD 和 GOS 均为取向梯度或晶内取向离散代理量。

织构采用 De la Vallee Poussin 核，半宽 5°，统一 5° 评估网格。像素加权 ODF 为主要状态描述，面积加权晶粒平均 ODF 用于权重敏感性。clipped_grid_entropy 在将负 ODF 网格值截断为零并归一至单位均值后计算；原始—去噪 ODF 距离使用单位均值谐波 ODF 的 L2 距离。

c 轴作为无方向轴处理：先令 AD 分量为正，若 AD 分量在 1e-12 容差内为零，则依次令 TD/RD、ND 分量为正。由此得到的 c 轴—AD 夹角为 0–90° 锐角，不能将峰分裂、方位重分配或织构展宽简化为单一净旋转。

轴向取向倾向以 AD 单轴拉伸、各滑移/孪生族最大绝对 Schmid 因子为几何假设；涵盖 basal <a>、prismatic <a>、pyramidal <a>、pyramidal <c+a> 及已验证的拉伸/压缩孪生族。未给定 CRSS，因而不据此判定实际激活顺序；Taylor 因子未计算，输出中的该列保持 NaN。

## 拉伸数据纳入规则

拉伸整合保留每条曲线原始 yield_status，并依据用户排除表确定是否纳入状态汇总。Rp0.2 仅汇总 yield_status=ok、数值有限且被用户接受的曲线。6.02 mm 状态的用户排除表仅接受 Y-1；Y-2 的 Rp0.2 标记为 not_reliable_no_elastic_segment_or_offset_crossing，虽其他拉伸字段仍为有限值，本流程仍依据用户排除表不将 Y-2 纳入状态汇总。因此该状态各汇总量的 n_valid=1，标准差为 0 仅表示单值汇总，不能解释为无离散性。
论文面向的拉伸汇总与 EBSD 相关性仅纳入 Rp0.2_MPa、UTS_engineering_MPa 和 strain_at_UTS_uniform_elongation_percent。max_true_stress_MPa 与 true_strain_at_max_true_stress_percent 由颈缩后均匀变形公式继续换算，不能作为局部真实应力—应变，故不纳入论文面向的汇总与相关性分析。

## EBSD—拉伸合并表的行粒度

`08_tensile_integration/ebsd_tensile_merged.csv` 是 EBSD 指标 × 拉伸曲线 × 拉伸量的笛卡尔长表，其自然键为 `(sample, ebsd_metric, tensile_sample, tensile_repeat, tensile_metric)`。同一 EBSD 值会随不同拉伸曲线和拉伸量重复，同一 tensile_value 以及 tensile_mean、tensile_sd、tensile_n_valid 也会随不同 EBSD 指标重复；这些重复值用于保留每个指标组合的可追溯性，不能作为独立实验观测或额外重复数。

`included_in_aggregate` 仅表示该条拉伸曲线的该项数值是否进入相应条件汇总；条件级均值、标准差和有效重复数以 `tensile_condition_aggregates.csv` 为准。相关性始终使用六个条件级状态，不按 `ebsd_tensile_merged.csv` 的行数扩增样本量。

## 结果的证据属性

| 属性 | 本流程中的结果 | 解释边界 |
| --- | --- | --- |
| 直接测量 | CTF 中的相、取向与质量字段，以及拉伸曲线得到的应力和应变指标 | 单个扫描区域或有限拉伸重复所代表的实验结果 |
| 取向派生量 | 晶粒形貌、晶界角度、织构和 c 轴分布 | 由注册重构与统计参数计算，不等同于变形机制的直接观测 |
| 代理指标 | KAM、GROD、GOS 及低角度晶界密度 | 反映取向梯度或亚结构特征，不是总位错密度的直接测量 |
| 推断 | 位错储存、动态回复、滑移或孪生对加工硬化的可能贡献 | 必须同时结合力学趋势、显微组织及文献；EBSD 与 Schmid/Taylor 结果不能唯一识别实际激活机制 |

## 原始—去噪比较与拉伸相关性

`07_raw_denoised_comparison` 给出每个注册汇总量的配对差值、相对差值、秩次一致性和相邻变形状态趋势一致性。若去噪改变相邻趋势方向或使带符号量跨越零值，结果会被标记为可能改变科学解释；原始数值仍保持主要地位。

`08_tensile_integration` 中的相关性采用六个直径状态的汇总值计算描述性 Spearman 系数，并报告逐一剔除状态后的范围。像素数未被视为实验重复，不进行高参数拟合；这些相关性用于提出假设，不能作为因果关系的证明。

## 拉伸相关性指标白名单

原始—去噪敏感性表保留全部数值型汇总量；拉伸相关性只纳入下列预先注册的科学指标，排除 grain_count、valid_pixel_count、total_area 和总长度等扫描支持量，并避免将同一 GROD/GOS 在不同 KAM 参数行中重复计算。

- 晶粒形貌：面积加权中值等效直径、面积加权中值长宽比、数目中值形状因子、面积加权长轴—AD 中值夹角。
- 晶界：晶界线密度、2–15° 低角度晶界长度分数、候选孪晶界长度分数。
- 晶内取向：主要参数（一阶邻域、5° 排除阈值）下的平均 KAM、平均 GROD 和面积加权平均 GOS。
- 织构：像素加权 texture index、M-index、c 轴—AD 平均夹角、c 轴位于 AD 30° 内比例及绕 AD 方位合成量。
- 轴向取向倾向：各滑移/孪生族面积加权平均最大绝对 Schmid 因子。

## 完整输出清单

- `README.md`
- `analysis_manifest.csv`

### `00_audit`
- `00_audit/scan_inventory.csv`
- `00_audit/raw_denoised_pair_audit.csv`
- `00_audit/raw_denoised_change_maps.png`

### `01_standard_maps`
- `01_standard_maps/7d_raw_maps.png`
- `01_standard_maps/7d_denoised_maps.png`
- `01_standard_maps/6.48d_raw_maps.png`
- `01_standard_maps/6.48d_denoised_maps.png`
- `01_standard_maps/6.02d_raw_maps.png`
- `01_standard_maps/6.02d_denoised_maps.png`
- `01_standard_maps/5.6d_raw_maps.png`
- `01_standard_maps/5.6d_denoised_maps.png`
- `01_standard_maps/5.25d_raw_maps.png`
- `01_standard_maps/5.25d_denoised_maps.png`
- `01_standard_maps/5d_raw_maps.png`
- `01_standard_maps/5d_denoised_maps.png`

### `02_grain_morphology`
- `02_grain_morphology/grain_morphology_by_grain.csv`
- `02_grain_morphology/grain_morphology_summary.csv`
- `02_grain_morphology/grain_morphology_trends.png`
- `02_grain_morphology/grain_morphology_quantiles.csv`（附加诊断）

### `03_boundaries`
- `03_boundaries/boundary_segments.csv`
- `03_boundaries/boundary_summary.csv`
- `03_boundaries/boundary_angle_distributions.csv`
- `03_boundaries/boundary_trends.png`
- `03_boundaries/boundary_detection_sensitivity.csv`（附加诊断）

### `04_intragranular`
- `04_intragranular/intragranular_by_pixel.csv`
- `04_intragranular/intragranular_by_grain.csv`
- `04_intragranular/intragranular_summary.csv`
- `04_intragranular/intragranular_trends.png`

### `05_texture`
- `05_texture/texture_summary.csv`
- `05_texture/c_axis_orientation_distribution.csv`
- `05_texture/pole_figures.png`
- `05_texture/inverse_pole_figures.png`
- `05_texture/odf_sections.png`
- `05_texture/texture_trends.png`
- `05_texture/c_axis_pole_figures_raw.png`
- `05_texture/c_axis_pole_figures_denoised.png`
- `05_texture/c_axis_mean_orientation_raw.png`
- `05_texture/c_axis_mean_orientation_denoised.png`
- `05_texture/texture_numerical_diagnostics.csv`
- `05_texture/raw_denoised_odf_distance.csv`
- `05_texture/task4_analysis_parameters.csv`
- `05_texture/texture_plot_color_limits.csv`

The raw and denoised c-axis figures are separate processing variants. The two pole-figure montages use one shared MRD scale, and the two mean-orientation trend figures use identical axis ranges.

### `06_axial_propensity`
- `06_axial_propensity/axial_propensity_by_grain.csv`
- `06_axial_propensity/axial_propensity_summary.csv`
- `06_axial_propensity/axial_propensity_trends.png`

### `07_raw_denoised_comparison`
- `07_raw_denoised_comparison/raw_denoised_metric_comparison.csv`
- `07_raw_denoised_comparison/raw_denoised_rank_agreement.csv`
- `07_raw_denoised_comparison/raw_denoised_metric_comparison.png`

### `08_tensile_integration`
- `08_tensile_integration/ebsd_tensile_merged.csv`
- `08_tensile_integration/ebsd_tensile_rank_correlations.csv`
- `08_tensile_integration/tensile_condition_aggregates.csv`
- `08_tensile_integration/ebsd_tensile_comparison.png`

### `09_sensitivity`（raw-only附加诊断）
- `09_sensitivity/morphology_sensitivity_summary.csv`
- `09_sensitivity/morphology_sensitivity_trends.png`

该附加模块使用 2°取向域的`minPixel=1/3/5/10`和15° HAGB晶粒的`minPixel=5`进行无平滑形貌对照；横向均为 AD。它用于区分低角度取向域细分与HAGB物理晶粒形貌，不把小域筛选差异解释为独立实验重复。

以上清单包含全部注册产物及生成器写出的附加诊断表。所有注册输入的 SHA-256、MATLAB/MTEX 版本和主要参数记录于 `analysis_manifest.csv`。
