# 图 4 的独立几何归一化版本

本目录为派生结果，不覆盖原图、原始 EBSD 数据或已有逐晶粒统计，也不替换稿件中的图片链接。

## 数据处理

- 来源：`results/mtex_grain_size_distribution_matrix/hagb_grain_size_by_grain.csv`，与原图 4 使用的统计汇总及直方图逐状态核对。
- 保持原始 EBSD 重构与筛选：HAGB 阈值 15°，晶粒至少含 5 个像素；保留原有接触视场边界的晶粒，不新增排除项。
- 按每个晶粒计算 `ECD_norm = ECD_observed × sqrt(d / 7.00)`，其中 d 为当前稿件的名义棒材直径，单位为 mm。
- 原始扫描 `6.48d` 对应当前稿件的名义 6.50 mm 状态，使用 6.50 mm 及由其计算的 13.78% 截面减缩率。派生表同时保留来源直径与名义直径，以明确对应关系。
- 归一化后重新分箱，六状态统一采用 2 μm 箱宽及原图的横轴范围；纵轴为晶粒数量频率，各状态总频率为 100%。柱高可能因重新分箱发生变化，不是简单移动原直方图。
- 黑线沿用原图的描述性对数正态拟合口径：对逐晶粒 ECD 的对数计算数量加权均值与总体标准差，再以 `100 × 箱宽 × PDF` 绘制。归一化后对数均值平移 `ln(sqrt(d/7.00))`，对数标准差不变。拟合不用于证明分布必然符合对数正态模型。
- 图中仅标注数量平均 ECD，不标注中位数。LAGB 和 2° 取向域未作归一化处理。

该换算基于均匀、轴对称、体积不变以及晶粒随宏观材料变形的假设，所得尺寸为换算至初始几何状态的纵截面比较指标，不是实测三维晶粒尺寸。

## 输出

- `hagb_grain_size_geometrically_normalized_six_state.png`：600 dpi 六状态组合图。
- 同名 `.tif`：600 dpi 无损压缩图。
- 同名 `.svg`：矢量图。
- `normalized_hagb_by_grain.csv`：逐晶粒归一化数据。
- `normalized_hagb_histograms.csv`：重新分箱后的数量频率。
- `normalized_hagb_summary.csv`：均值、归一化因子及拟合参数。
- `provenance.json`：来源校验值、映射关系与验证记录。

## 建议图注

不同旋锻状态下 HAGB 晶粒几何归一化等效圆直径的数量频率分布。以 15° 取向差阈值重构晶粒，并保留至少 5 个像素的晶粒；按 ECD_norm = ECD_observed × sqrt(d/d₀) 将纵截面晶粒尺寸换算至初始几何状态，其中 d₀ = 7.00 mm，d 为各状态名义棒材直径。各状态采用相同的 2 μm 分箱宽度，黑线为描述性对数正态拟合，图中标注值为归一化后的数量平均 ECD。

## 复现

运行脚本：`tools/plot_normalized_hagb_grain_size.py`。

本次在独立临时依赖目录安装 Matplotlib，不修改项目或系统 Python 环境：

```bash
PYTHONPATH=.codex_tmp/hagb_plot_python_deps /home/abcd/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 tools/plot_normalized_hagb_grain_size.py
```
