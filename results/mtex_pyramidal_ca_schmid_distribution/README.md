# 锥面 `<c+a>` 滑移施密特因子分布

本目录给出六种冷变形状态在 AD 单轴拉伸假设下，锥面 `<c+a>` 滑移族最大绝对施密特因子的晶粒面积加权分布。每个晶粒均对 MTEX `slipSystem.pyramidalCA` 的全部对称等价变体取最大绝对值；raw CTF 为主结果，denoised CTF 用于敏感性检查。

## 解释边界

施密特因子仅描述给定加载方向下的几何分切条件。该统计未引入 CRSS、旋锻多轴应力、邻晶约束、局部应力和滑移迹线，因此不能单独用于证明锥面 `<c+a>` 已经启动。结合项目已有 CRSS 归一化结果，更稳妥的表述是：随冷变形增加，锥面 `<c+a>` 的几何有利程度相对增强，但柱面 `<a>` 仍具有更低的名义启动阻力。

## 文件

- `pyramidal_ca_schmid_distribution_6state.png/.tif/.pdf`：六状态分布图。
- `pyramidal_ca_schmid_distribution.csv`：逐分箱面积加权频率。
- `pyramidal_ca_schmid_summary.csv`：各状态均值、中位数和高因子面积比例。
- `pyramidal_ca_schmid_parameters.csv`：计算参数与解释边界。

复现命令：

```bash
tools/run_mtex_guarded.sh --timeout 3600 \
  "generate_pyramidal_ca_schmid_distribution(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')),string(fullfile(pwd,'results','mtex_pyramidal_ca_schmid_distribution')))"
```
