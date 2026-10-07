# MTEX manuscript evidence

该目录保存当前正文表 3 所依赖的稳定 Schmid 因子汇总和趋势图。结果于 2026-08-24 在 Linux MATLAB R2024a、MTEX 6.1.1 环境中，由受控启动器生成：

```bash
tools/run_mtex_guarded.sh --timeout 1800 \
  "generate_comprehensive_axial_propensity(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')),string(fullfile(pwd,'results','mtex_manuscript_evidence')))"
```

## 文件

- `06_axial_propensity/axial_propensity_summary.csv`：六状态 raw/denoised、四类滑移族和两类孪生族的条件级汇总。正文以 raw 行为主。
- `06_axial_propensity/axial_propensity_trends.png`：raw/denoised 趋势图。
- `axial_propensity_by_grain.csv`：生成器会创建的逐晶粒中间表，约 65 MB，可由 CTF 重建，故不纳入 Git。

计算假设为沿 AD 的单轴拉伸几何分切。未输入 CRSS；孪生只计算最大绝对 Schmid 因子，不能解析极性。因此这些结果不能单独判定实际滑移或孪生启动。

本次受控批处理的 MATLAB 退出码为 0，未触发超时，结束后运行 SID 无残留进程。
