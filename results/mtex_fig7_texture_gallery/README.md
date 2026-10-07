# Fig. 7 参考式织构图

本目录给出 Gr4B23271 商业纯钛六种冷变形状态的极图（PF）和逆极图（IPF）。输入为 `data/ebsd_kpl_250221_7_df/scans/` 下六组原始 CTF 的 Ti-Hex 相；未使用去噪 CTF 替代定量主数据。

## 输出

- `reference_style/01_state_figures/`：六张独立状态小图，按 0、14.31%、26.04%、36.00%、43.75% 和 48.98% 冷变形量排列。
- `reference_style/montages/fig7_reference_style_six_state_montage.{png,tif,pdf}`：六状态组合大图。
- `01_state_figures/` 和 `montages/fig7_texture_six_state_montage.*`：原有统一色标版本，用于跨状态按颜色直接比较。
- `fig7_texture_summary.csv`：输入路径、有效取向数、各状态自动色标上限、统一色标上限及两套小图路径。
- `fig7_texture_parameters.csv`：ODF 与绘图参数。
- `fig7_reference_style_qa.md`：与参考 Fig. 7 的逐项形式对照及无等高线验证记录。

## 图形与参数约定

- PF：`{0001}`、`{10-10}`、`{11-20}`；试样坐标为 `AD`、`TD/RD`，圆心为 `ND`。
- IPF：`AD`、`TD/RD`、`ND` 三个方向。
- ODF 核半宽 5°；统一色标版绘图采样分辨率为 5°，参考式平滑版为 2°；输出 600 dpi。
- 两套正式输出均采用 `smooth` 连续密度填色，不绘制等高线；参考式 2°采样进一步降低分级带状感。
- 参考式小图采用各状态、各行独立 MRD 上限，并显式标注 `Max` 与 `Min=0.00`；它适合显示峰位和分布形态，但不能仅按颜色比较不同状态的绝对强度。
- 统一色标版的 PF/IPF 上限均为 7 MRD，用于跨状态比较。

版式参考 Xia 等发表于 *Materials Science & Engineering A* 975 (2026) 150869 的 Fig. 7，仅复用 PF/IPF 的内容层级与构图逻辑；所有密度场均由本项目 EBSD 数据重新计算。生成器为 `tools/mtex/generate_fig7_texture_gallery.m`。
