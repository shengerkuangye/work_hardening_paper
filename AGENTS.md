# Gr4B23271 Work Hardening Paper Agent

## Project Mission

This project supports writing an academic paper focused on the work hardening mechanism of Gr4B23271 commercially pure titanium during cold deformation. Treat the repository as a research workspace: preserve raw data, cite evidence carefully, and separate measured results from interpretation.

## Core Research Focus

- Main material: Gr4B23271 commercially pure titanium.
- Main topic: cold-deformation-induced work hardening mechanism.
- Analytical sequence: diameter reduction -> cold deformation ratio -> tensile response -> work hardening behavior -> microstructure/EBSD evidence -> mechanism discussion.
- Preferred mechanism language: dislocation accumulation, grain orientation/texture evolution, grain boundary effects, possible twinning or slip-system activity only when supported by EBSD or literature evidence.

## Key Local Sources

- Tensile summaries:
  - `data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/gr4b23271_cold_deformation_tensile_summary.csv`
  - `data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/gr4b23271_tensile_final_by_diameter_user_exclude.csv`
  - `data/tensile_data/gr4b23271_lab_tensile_by_diameter.csv`
- Raw tensile data:
  - `data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/`
- Microstructure and EBSD outputs:
  - `data/metallography_2025_02_18/images/`
  - `data/ebsd_kpl_250221_7_df/scans/`
- Extracted literature notes:
  - `ti6554_cryogenic_mechanical_properties_ml_extracted.txt`
  - `ml_constitutive_modelling_material_nonlinearity_review_extracted.txt`
  - `arrhenius_bp_hot_deformation_hypereutectoid_steel_extracted.txt`
- Original references:
  - `references/`

## Working Rules

- Do not overwrite raw experimental data. Create derived tables, figures, and notes as new files with clear names.
- When analyzing mechanical data, document exclusions and reliability flags. In particular, preserve distinctions such as `ok` versus `not_reliable_no_elastic_segment_or_offset_crossing`.
- Use SI units consistently. Report stress in MPa, strain as dimensionless or percent as appropriate, and deformation ratio in percent.
- For paper writing, prefer precise academic Chinese or English. Avoid unsupported claims such as "obviously", "proves", or "the only mechanism".
- When discussing mechanisms, link each claim to one of:
  - tensile trend,
  - work hardening rate curve,
  - EBSD/map evidence,
  - optical microstructure,
  - cited literature.
- If a conclusion needs EBSD-derived quantities such as KAM, GND density, grain size, misorientation, or texture intensity, check whether those data exist locally before asserting them.
- On Windows, run MATLAB batches hidden with `Start-Process -WindowStyle Hidden` and `-nodesktop -nosplash -noFigureWindows -batch`; after exit, inventory `matlab.exe`, `MATLAB.exe`, and `MATLABWindow.exe`, clean only exact run-owned orphan `MATLABWindow.exe` PIDs after verifying launcher/main counts are zero, and never kill a pre-existing or interactive user process.

## Academic Wording Requirements

- Manuscript text must use formal academic language, not planning-note or agent-workflow language.
- Do not write phrases such as "核心证据链", "证据链如下", "本文的故事线", "主线是", "异常点", "打通逻辑", "说明一下", or similar informal/project-management wording in the paper body.
- Prefer manuscript-ready phrasing such as:
  - "本文围绕旋锻变形量、拉伸响应与显微组织演化之间的关系展开研究。"
  - "为阐明旋锻冷变形对 Gr4B23271 商业纯钛强塑性匹配的影响，本文系统分析了不同变形量下的力学性能和显微组织特征。"
  - "结合拉伸曲线、加工硬化行为、金相组织和 EBSD 表征结果，对强度提升及塑性变化的组织机制进行讨论。"
  - "对于偏离整体趋势的实验结果，应结合原始曲线可靠性、金相组织和 EBSD 指标进行分析。"
- Distinguish drafting notes from manuscript prose. Terms such as "evidence chain", "to-do", "data gap", "check EBSD", and "possible explanation" may appear in notes or outlines, but should be rewritten before being used in Introduction, Results, Discussion, or Conclusions.
- Use cautious causal language unless the evidence is direct: prefer "表明", "说明", "可能与...有关", "可归因于", "主要受...影响"; avoid "证明", "必然导致", "唯一原因", "显然".
- For abnormal or non-monotonic data, use academic wording such as "偏离整体变化趋势的结果", "非单调变化", "局部差异", or "离散性增加" rather than casual terms such as "异常点" in final manuscript text.

## Suggested Paper Structure

1. Introduction: Gr4B23271 applications, cold deformation strengthening, current gap in connecting mechanical response with microstructure.
2. Experimental Materials and Methods: material, cold drawing/deformation route, tensile testing, metallography, EBSD.
3. Results:
   - deformation ratio and sample grouping,
   - tensile properties versus deformation ratio,
   - true stress-strain and work hardening response,
   - microstructure and EBSD observations.
4. Discussion:
   - work hardening stages,
   - dislocation storage and dynamic recovery balance,
   - texture/grain-boundary contribution,
   - strength-ductility tradeoff.
5. Conclusions: concise, data-backed points.

## Agent Behavior

When asked to help with this project:

- First inspect the relevant local data or text instead of relying on memory.
- Keep edits scoped to the requested manuscript, analysis, figure, or note.
- Prefer reproducible analysis scripts or clearly named derived CSV files for numerical work.
- In final responses, summarize what changed and mention any remaining uncertainty or data gap.
