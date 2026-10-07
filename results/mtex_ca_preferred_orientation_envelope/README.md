# 锥面 `<c+a>` 优先启动的c轴—AD取向包络

本目录使用 Xiong 等商业纯钛 Grade 4 的CRSS，计算在AD单轴拉伸下，什么样的c轴—AD夹角可能使锥面 `<c+a>` 比柱面和基面 `<a>` 具有更高的名义启动倾向。

## 判据与扫描参数

定义：

\[
R=\frac{m_{\langle c+a\rangle}/756}
{\max(m_{\mathrm{prismatic}}/154,\;m_{\mathrm{basal}}/252)}.
\]

`R > 1` 表示锥面 `<c+a>` 的CRSS归一化驱动力高于柱面和基面 `<a>`。每个滑移族均对全部对称等价变体取最大绝对施密特因子。

- c轴—AD锐角：0°–25°，步长0.1°。
- 绕c轴方位角：0°–359.75°，步长0.25°。
- 取向晶格：251 × 1440，共361440个加载方向。
- CRSS：柱面 `<a>` 154 MPa、基面 `<a>` 252 MPa、锥面 `<c+a>` 756 MPa。

## 结果

| c轴—AD夹角 | `<c+a>` 优先情况 |
|---:|---|
| 0°–9.1° | 所有绕c轴方位均满足 `R > 1` |
| 9.2°–10.6° | 仅部分绕c轴方位满足 `R > 1` |
| ≥10.7° | 没有绕c轴方位满足 `R > 1` |

部分方位区间内，满足 `<c+a>` 优先的方位比例迅速降低：9.2°为80.4%，9.5°为50.4%，10.0°为23.8%，10.5°为4.6%，10.6°仅为1.25%。由于角度扫描步长为0.1°，连续边界应表述为“全方位边界位于9.1°–9.2°之间，可能方位边界位于10.6°–10.7°之间”，不宜将离散端点作为更高精度的材料常数。

因此，若只利用c轴—AD夹角筛选：

- `≤9.1°` 可以判为当前CRSS假设下的稳健 `<c+a>` 优先取向；
- `9.2°–10.6°` 必须进一步考虑绕c轴方位，不能只凭夹角判断；
- `≥10.7°` 在当前单轴和固定CRSS模型中不能优先于已纳入的两类 `<a>` 滑移。

## 与项目EBSD的对照

本项目六种raw状态的保留晶粒平均取向中，最小c轴—AD夹角为21.20°–36.70°，显著超出10.6°的可能优先边界。因此，当前视场中 `<c+a>` 优先晶粒为0与取向包络计算一致。

需要注意，`<c+a>` 施密特因子本身可在更大夹角处达到接近0.5，但这不等于优先启动。随着c轴偏离AD，低CRSS的基面或柱面 `<a>` 获得足够的分切驱动力，`R` 会降至1以下。

## 适用边界

该结果是以参考文献CRSS为固定输入的取向空间筛选，不是完整晶体塑性预测。局部多轴应力、邻晶约束、残余应力、柱面滑移加工硬化和各CRSS随变形历史的变化均可能扩展或收缩实际 `<c+a>` 活动区间。由于参考文献没有提供独立的锥面 `<a>` CRSS，本计算未将锥面 `<a>` 纳入竞争集合。

## 文件

- `ca_preferred_orientation_envelope.csv`：逐0.1°夹角的优先方位比例和驱动力比范围。
- `ca_preferred_orientation_parameters.csv`：CRSS、扫描分辨率和两类边界。
- `ca_preferred_orientation_envelope.png`、`.tif`：取向包络图。
- 生成器：`tools/mtex/generate_ca_preferred_orientation_envelope.m`。

复现命令：

```bash
tools/run_mtex_guarded.sh --timeout 1200 \
  "scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); outputRoot=string(fullfile(pwd,'results','mtex_ca_preferred_orientation_envelope')); generate_ca_preferred_orientation_envelope(scanRoot,outputRoot)"
```
