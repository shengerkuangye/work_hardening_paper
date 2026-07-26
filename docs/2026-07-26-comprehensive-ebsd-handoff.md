# Gr4B23271 综合 EBSD 分析交接

## 1. 分支与范围

- 工作分支：`feature/comprehensive-ebsd`
- 基线分支：`main`
- 远端目标：`origin/feature/comprehensive-ebsd`
- 试样坐标：图像横向 `x = AD`（棒材轴向），`y = TD/RD`，表面法向 `z = ND`
- 状态：7、6.48、6.02、5.6、5.25 和 5 mm，对应冷变形量约 0、14.31%、26.04%、36.00%、43.75% 和 48.98%
- 定量主结果采用 raw 数据；denoised 数据用于处理敏感性检查。

该分支包含 EBSD 原始/降噪审计、标准图、取向域与晶粒形貌、晶界错配、KAM/GROD/GOS、织构、c 轴取向分布、轴向加载几何倾向、拉伸数据联合表以及阈值敏感性分析。

## 2. 主要科学结论

1. 六个状态的像素加权 c 轴—AD 平均夹角约为 78.6°–79.8°，97.1%–99.7% 的 c 轴与 AD 夹角大于 60°。现有数据不支持 c 轴整体持续向 AD 单调旋转。
2. 完整 `{0001}` c 轴极图显示主导织构组分随变形发生非单调迁移和竞争。更稳妥的表述是“织构组分重新分配和取向选择”，而不是把离散试样之间的主峰变化解释为同一批晶粒连续旋转。
3. KAM、GROD、GOS 和 LAGB 在首个变形状态已经明显增加。因此，“后期才发生位错增殖”的严格先后模型不符合当前 EBSD 观察。
4. 36%–43.75% 区间的低角度边界密度、取向域细分和轴向排列增强，更适合解释为已有取向梯度向低角度变形边界网络的组织化。
5. 2°取向域的细化不能等同于 HAGB 晶粒细化或再结晶；15° HAGB 晶粒主要表现为沿 AD 拉长，而没有呈现连续晶粒尺寸降低。
6. raw、denoised 和共同测点支持区域给出一致的 c 轴结论，但每个状态只有一个 EBSD 视场，非单调局部变化不能直接作为机制转折或因果证据。

推荐的论文表述为：

> 旋锻初期即伴随晶内取向梯度和低角度变形边界形成；随着变形量提高，取向梯度进一步组织为更细密且轴向排列的低角度取向域，同时织构组分发生非单调的重新分配。现有 EBSD 结果未显示 c 轴整体向棒材轴向持续旋转。

## 3. 关键源码与入口

- 综合入口：`tools/mtex/run_comprehensive_ebsd_analysis.m`
- 输出合同：`tools/mtex/comprehensive_ebsd_output_contract.m`
- 晶内与织构：`tools/mtex/generate_comprehensive_intragranular_texture.m`
- c 轴一维和球面分布：`tools/mtex/generate_c_axis_distribution_functions.m`
- 形貌敏感性：`tools/mtex/generate_comprehensive_morphology_sensitivity.m`
- raw/denoised 对比：`tools/mtex/generate_comprehensive_raw_denoised_comparison.m`
- 拉伸联合分析：`tools/mtex/generate_comprehensive_tensile_integration.m`
- 机制与补实验路线图：`docs/2026-07-23-ebsd-mechanism-and-paper-roadmap.md`

## 4. 运行方法

依赖 MATLAB R2025a 和 MTEX 6.1.1。进入项目根目录后运行：

```matlab
run('C:/Users/22069/Documents/MATLAB/mtex-6.1.1/startup_mtex.m');
addpath('tools/mtex');
projectRoot = string(pwd);
outputRoot = fullfile(projectRoot,'results','mtex_ebsd_comprehensive');
run_comprehensive_ebsd_analysis(projectRoot,outputRoot);
```

只在已有结果上重写清单和说明并复核输入哈希及输出合同：

```matlab
run_comprehensive_ebsd_analysis(projectRoot,outputRoot, ...
  struct('finalize_only',true));
```

## 5. 本地产物

完整产物位于：

`results/mtex_ebsd_comprehensive/`

其中 c 轴定量表和图位于：

`results/mtex_ebsd_comprehensive/05_texture/`

重点文件包括：

- `c_axis_ad_distribution_function.csv`
- `c_axis_spherical_distribution_function.csv`
- `c_axis_distribution_parameters.csv`
- `c_axis_ad_distribution_raw.png`
- `c_axis_ad_distribution_denoised.png`
- `c_axis_spherical_distribution_raw.png`
- `c_axis_spherical_distribution_denoised.png`

完整结果约 2.14 GB，且包含多个超过 GitHub 100 MB 限制的派生 CSV，因此不纳入 Git 提交。所有原始 CTF 和拉伸输入均保持不变，结果可通过上述入口重新生成。

## 6. 验证

交接前应至少运行：

```matlab
test_comprehensive_ebsd_contract;
test_comprehensive_ebsd_analysis;
test_comprehensive_intragranular_texture;
test_comprehensive_morphology_sensitivity;
test_c_axis_distribution_functions;
```

c 轴真实产物还需满足：

- 每组一维分布概率和为 1；
- 球面分布 `sum(cell_weight .* MRD) = 1`；
- MRD 有限且非负；
- raw 和 denoised 各覆盖六个状态；
- 四张图的坐标、色标和试样方向标注一致。

## 7. 尚需补充的实验

若目标是建立高可信度加工硬化机制，优先补充：

1. 每个变形量至少 3 个独立 EBSD 视场，并覆盖径向和周向位置，用于区分真实非单调演化与抽样差异。
2. TEM/STEM 位错与亚晶结构，用于直接验证位错胞、位错墙和低角度亚边界。
3. XRD 线宽或 HR-EBSD，用于获得相对或绝对位错密度及晶格应变。
4. 孪晶界识别、滑移痕迹及 Schmid/CRSS 联合分析，用于约束实际变形模式；仅凭 c 轴极图不能判定主导滑移系。
5. 重新核算均匀塑性区间内的真应力—真塑性应变和加工硬化率，并报告重复试样离散性。

## 8. 工作区注意事项

本 worktree 仍保留以下未纳入本次 EBSD 交接提交的独立改动：

- `tools/plot_manuscript_figures.py`
- `tools/test_plot_manuscript_figures.py`

它们涉及拉伸真塑性应变和加工硬化图的修订，应在独立测试和审查后单独提交，避免与 EBSD 分支收尾混合。
