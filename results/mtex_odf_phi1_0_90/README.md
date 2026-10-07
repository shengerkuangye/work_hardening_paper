# 图7 ODF横坐标0–90°版本

本目录将原图的 φ₁ 显示区间从0–360°改为0–90°。保留原始Ti-Hex取向、试样对称性1、晶体对称性6/mmm、5°核半宽和统一MRD色标。未进行象限折叠、附加试样对称化或区间重新归一化。原0–360°图和原始CTF数据均保留。

各子图横轴φ₁与纵轴Φ均为0–90°，绘图区宽高比及数据单位比例均为1∶1。脚本逐一校验18个绘图区的像素宽高相等，Word插图按原图纵横比等比例缩放。

- `odf_selected_phi2_phi1_0_90.png`：六状态×三个φ₂截面的绘图输出。
- `odf_selected_phi2_phi1_0_90.pdf`：同图矢量输出。
- `density_*_phi2_*.csv`：1°网格的原始计算密度，行对应Φ=0:90°，列对应φ₁=0:90°。
- `odf_window_peak_positions.csv`：仅0–90°显示区间内的截面峰值及坐标。
- `odf_full_space_summary.csv`：重新计算的全取向空间统计。
- `odf_full_section_peaks.csv`：0–359°完整截面统计，便于核对区间变化。
- `odf_models.mat`：本次计算的ODF对象。
- `plot_provenance.json`：色标、采样间隔、绘图处理及密度文件校验值。

图中Max ODF仍表示全取向空间最大值，不表示缩小后显示窗口的最大值。绘图取5°采样，与原图分辨率一致；负密度仅在显示时截为0，与原绘图处理一致，CSV保留计算值。

计算入口：`tools/mtex/export_odf_phi1_0_90.m`；必须通过`tools/run_mtex_guarded.sh`运行。绘图入口：`tools/plot_odf_phi1_0_90.py`，依赖NumPy、pandas和Matplotlib。当前Word图7及对应段落已按窗口内峰值调整，其他图像保持不变。
