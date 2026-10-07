# Fig. 14参考形式：六直径KAM-GND图

## 图件内容

- 三行分别比较7-6.48 mm、6.02-5.6 mm和5.25-5 mm。
- 每行左侧为两张GND空间图，右侧为对应的上下两组像素频率分布。
- 六张空间图采用统一色标，不绘制等高线。

## 计算方法

- 输入为raw CTF；KAM采用一阶邻域和5 deg排除阈值。
- GND按rho_GND=2*KAM/(step*b)逐像素换算，其中KAM使用弧度。
- EBSD步长为0.5 um，Burgers矢量取0.295 nm。
- 频率分布按有效Ti-Hex像素计数。

## 解释边界

该结果是基于KAM的几何必要位错密度估算，反映局部取向梯度，不等同于包含统计存储位错在内的总位错密度。结果受步长、邻域阶数、阈值、索引质量和数据清理方法影响。

## 建议图注

不同冷变形状态下Gr4B23271商业纯钛的KAM估算GND密度空间分布及统计结果：（a-c）7和6.48 mm样品；（d-f）6.02和5.6 mm样品；（g-i）5.25和5 mm样品。GND密度由一阶邻域KAM按rho_GND=2*KAM/(step*b)估算，EBSD步长为0.5 um，b=0.295 nm。空间图采用统一色标，统计分布按有效Ti-Hex像素计数。

## 输出

- `reference_style_fig14_gnd_six_diameters.*`：600 dpi PNG、TIFF和PDF。
- `reference_style_fig14_gnd_summary.csv`：六状态KAM和GND统计。
- `reference_style_fig14_gnd_distribution.csv`：完整分布数据。
- `reference_style_fig14_gnd_parameters.csv`：计算参数和解释边界。
