function run_comprehensive_ebsd_analysis(projectRoot, outputRoot, options)
%RUN_COMPREHENSIVE_EBSD_ANALYSIS Generate the registered derived bundle.

arguments
  projectRoot (1,1) string
  outputRoot (1,1) string
  options (1,1) struct = struct()
end

assert(isfolder(projectRoot), "Project root not found: %s", projectRoot);
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
projectRoot = canonical_path(projectRoot);
outputRoot = validate_output_root(projectRoot, outputRoot);
scanRoot = fullfile(projectRoot, "data", "ebsd_kpl_250221_7_df", ...
  "scans");
assert(isfolder(scanRoot), "Registered EBSD scan root is missing.");
catalog = comprehensive_ebsd_catalog(scanRoot);
hashesBefore = hash_catalog_inputs(catalog);
tensileInputs = registered_tensile_inputs(projectRoot);
tensileHashesBefore = hash_paths(tensileInputs.absolute_path);
finalizeOnly = isfield(options, "finalize_only") && options.finalize_only;

% outputRoot has been resolved, constrained to this project, and excluded
% from source-data/reference locations before this fresh-tree operation.
contract = comprehensive_ebsd_output_contract();
if finalizeOnly
  assert(isfolder(outputRoot), ...
    "finalize_only requires an existing derived-output directory.");
  assert_owned_output_root(outputRoot);
  assert(isfile(fullfile(outputRoot, "analysis_manifest.csv")), ...
    "finalize_only requires an existing provenance manifest.");
  verify_manifest_hashes(outputRoot, catalog, hashesBefore, ...
    tensileInputs, tensileHashesBefore);
else
  tensileOptions = load_registered_tensile_tables(tensileInputs);
  if isfolder(outputRoot)
    assert_owned_output_root(outputRoot);
    rmdir(outputRoot, "s");
  end
  mkdir(outputRoot);
  write_output_marker(outputRoot);
  for directory = contract.directories'
    mkdir(fullfile(outputRoot, directory));
  end

  generate_comprehensive_ebsd_audit(scanRoot, outputRoot);
  generate_comprehensive_maps_morphology_boundaries(scanRoot, outputRoot);
  generate_comprehensive_intragranular_texture(scanRoot, outputRoot);
  generate_c_axis_distribution_functions( ...
    fullfile(outputRoot, "05_texture"));
  generate_comprehensive_axial_propensity(scanRoot, outputRoot);
  [comparison, ~] = ...
    generate_comprehensive_raw_denoised_comparison(outputRoot);
  generate_comprehensive_tensile_integration(projectRoot, outputRoot, ...
    comparison, tensileOptions);
end

hashesAfter = hash_catalog_inputs(catalog);
assert(isequal(hashesBefore, hashesAfter), ...
  "One or more registered CTF inputs changed during analysis.");
tensileHashesAfter = hash_paths(tensileInputs.absolute_path);
assert(isequal(tensileHashesBefore, tensileHashesAfter), ...
  "One or more registered tensile inputs changed during analysis.");
if ~finalizeOnly
  write_manifest(outputRoot, catalog, hashesBefore, tensileInputs, ...
    tensileHashesBefore, contract);
end
verify_manifest_hashes(outputRoot, catalog, hashesBefore, ...
  tensileInputs, tensileHashesBefore);
write_manuscript_readme(outputRoot, contract, catalog);
verify_output_contract(outputRoot, contract, height(catalog));
fprintf("Comprehensive EBSD bundle written to %s\n", outputRoot);
end

function markerPath = output_marker_path(outputRoot)
markerPath = fullfile(outputRoot, ".comprehensive_ebsd_owned");
end

function write_output_marker(outputRoot)
markerPath = output_marker_path(outputRoot);
fileId = fopen(char(markerPath), "w", "n", "UTF-8");
assert(fileId >= 0, "Could not write output ownership marker.");
cleanupFile = onCleanup(@() fclose(fileId));
fprintf(fileId, "gr4b23271-comprehensive-ebsd-v1\n");
end

function assert_owned_output_root(outputRoot)
markerPath = output_marker_path(outputRoot);
assert(isfile(markerPath), ...
  "Refusing to reuse or delete an unowned existing outputRoot.");
marker = strtrim(string(fileread(markerPath)));
assert(marker == "gr4b23271-comprehensive-ebsd-v1", ...
  "Output ownership marker is invalid.");
end

function outputRoot = validate_output_root(projectRoot, outputRoot)
outputRoot = canonical_path(outputRoot);
projectPrefix = lower(projectRoot + string(filesep));
assert(startsWith(lower(outputRoot), projectPrefix), ...
  "outputRoot must be a derived-output folder inside projectRoot.");
assert(lower(outputRoot) ~= lower(projectRoot), ...
  "outputRoot cannot be projectRoot.");
allowed = false;
for allowedName = ["results",".codex_tmp"]
  allowedRoot = lower(canonical_path(fullfile(projectRoot, allowedName)));
  allowed = allowed || startsWith(lower(outputRoot), ...
    allowedRoot + string(filesep));
end
assert(allowed, ...
  "outputRoot must be below results/ or .codex_tmp/.");
for protectedName = ["data","references"]
  protectedPath = lower(canonical_path(fullfile(projectRoot, ...
    protectedName)));
  assert(lower(outputRoot) ~= protectedPath && ...
    ~startsWith(lower(outputRoot), protectedPath + string(filesep)), ...
    "outputRoot cannot be inside protected source folder %s.", ...
    protectedName);
end
assert(strlength(outputRoot) > strlength(projectRoot) + 2, ...
  "outputRoot is too broad for a fresh derived-output tree.");
end

function path = canonical_path(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end

function hashes = hash_catalog_inputs(catalog)
hashes = strings(height(catalog), 1);
for inputIndex = 1:height(catalog)
  hashes(inputIndex) = sha256_file(catalog.input_path(inputIndex));
end
end

function inputs = registered_tensile_inputs(projectRoot)
name = ["tensile_detail";"tensile_user_exclusion_summary"];
relative_path = [ ...
  "data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/" + ...
    "gr4b23271_cold_deformation_tensile_summary.csv"; ...
  "data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/" + ...
    "gr4b23271_tensile_final_by_diameter_user_exclude.csv"];
absolute_path = strings(2,1);
for inputIndex = 1:2
  absolute_path(inputIndex) = canonical_path(fullfile(projectRoot, ...
    replace(relative_path(inputIndex), "/", filesep)));
  assert(isfile(absolute_path(inputIndex)), ...
    "Registered tensile input is missing: %s", absolute_path(inputIndex));
end
inputs = table(name, relative_path, absolute_path);
end

function hashes = hash_paths(paths)
hashes = strings(numel(paths), 1);
for inputIndex = 1:numel(paths)
  hashes(inputIndex) = sha256_file(paths(inputIndex));
end
end

function tensileOptions = load_registered_tensile_tables(inputs)
tensileOptions = struct();
tensileOptions.tensile_detail = readtable(inputs.absolute_path(1), ...
  "TextType", "string", "VariableNamingRule", "preserve");
tensileOptions.user_summary = readtable(inputs.absolute_path(2), ...
  "TextType", "string", "VariableNamingRule", "preserve");
end

function write_manifest(outputRoot, catalog, hashes, tensileInputs, ...
  tensileHashes, contract)
category = strings(0,1);
name = strings(0,1);
value = strings(0,1);
path = strings(0,1);
sha256 = strings(0,1);
evidence_class = strings(0,1);
manifest = table(category, name, value, path, sha256, evidence_class);

for inputIndex = 1:height(catalog)
  row = table("input_ctf", ...
    catalog.sample(inputIndex) + "_" + catalog.variant(inputIndex), ...
    "registered native-grid input", catalog.input_path(inputIndex), ...
    hashes(inputIndex), ...
    input_evidence_class(catalog.variant(inputIndex)), ...
    'VariableNames', manifest.Properties.VariableNames);
  manifest = [manifest; row]; %#ok<AGROW>
end
for inputIndex = 1:height(tensileInputs)
  row = table("input_tensile", tensileInputs.name(inputIndex), ...
    "registered tensile integration input", ...
    tensileInputs.relative_path(inputIndex), tensileHashes(inputIndex), ...
    "derived tensile summary input with preserved reliability/exclusion fields", ...
    'VariableNames', manifest.Properties.VariableNames);
  manifest = [manifest; row]; %#ok<AGROW>
end
manifest = [manifest; table("software", "MATLAB", string(version), ...
  "", "", "analysis environment", 'VariableNames', ...
  manifest.Properties.VariableNames)];
manifest = [manifest; table("software", "MTEX", ...
  string(getMTEXpref("version")), "", "", "analysis environment", ...
  'VariableNames', manifest.Properties.VariableNames)];
scientificRows = [ ...
  "direct_measurement", "CTF_and_tensile_fields", ...
    "measured phase/orientation/quality fields and tensile stress-strain quantities", ...
    "direct measurement"; ...
  "orientation_derived", "morphology_boundary_texture", ...
    "quantities calculated from registered EBSD reconstruction and orientations", ...
    "orientation-derived result"; ...
  "proxy_indicator", "KAM_GROD_GOS_LAGB", ...
    "orientation-gradient and substructure proxies; not total dislocation density", ...
    "proxy indicator"; ...
  "inference_limit", "mechanism_interpretation", ...
    "mechanism attribution requires combined evidence and is not uniquely identified", ...
    "inference, not direct measurement"];
for scienceIndex = 1:size(scientificRows, 1)
  manifest = [manifest; table(scientificRows(scienceIndex, 1), ...
    scientificRows(scienceIndex, 2), scientificRows(scienceIndex, 3), ...
    "", "", scientificRows(scienceIndex, 4), 'VariableNames', ...
    manifest.Properties.VariableNames)]; %#ok<AGROW>
end

parameterNames = string(fieldnames(contract.parameters));
for parameterIndex = 1:numel(parameterNames)
  parameterName = parameterNames(parameterIndex);
  row = table("registered_parameter", parameterName, ...
    value_to_string(contract.parameters.(parameterName)), "", "", ...
    "analysis parameter", 'VariableNames', ...
    manifest.Properties.VariableNames);
  manifest = [manifest; row]; %#ok<AGROW>
end
writetable(manifest, fullfile(outputRoot, "analysis_manifest.csv"));
end

function result = input_evidence_class(variant)
if string(variant) == "raw"
  result = "direct measured input; primary quantitative source";
else
  result = "processed paired sensitivity input; not an independent replicate";
end
end

function result = value_to_string(value)
result = serialize_manifest_value(value);
end

function verify_manifest_hashes(outputRoot, catalog, hashes, ...
  tensileInputs, tensileHashes)
manifestPath = fullfile(outputRoot, "analysis_manifest.csv");
importOptions = delimitedTextImportOptions("NumVariables", 6);
importOptions.DataLines = [2 Inf];
importOptions.Delimiter = ",";
importOptions.VariableNames = ...
  {'category','name','value','path','sha256','evidence_class'};
importOptions.VariableTypes = repmat("string", 1, 6);
importOptions.ExtraColumnsRule = "error";
importOptions.EmptyLineRule = "read";
manifest = readtable(manifestPath, importOptions);
required = ["category","name","path","sha256"];
assert(all(ismember(required, ...
  string(manifest.Properties.VariableNames))), ...
  "Analysis manifest lacks registered provenance columns.");
inputRows = manifest.category == "input_ctf";
assert(nnz(inputRows) == height(catalog), ...
  "Analysis manifest must contain one row per registered CTF input.");
for inputIndex = 1:height(catalog)
  expectedName = catalog.sample(inputIndex) + "_" + ...
    catalog.variant(inputIndex);
  row = inputRows & manifest.name == expectedName;
  assert(nnz(row) == 1, ...
    "Analysis manifest input key is missing or duplicated: %s", ...
    expectedName);
  assert(manifest.path(row) == catalog.input_path(inputIndex), ...
    "Analysis manifest path mismatch for %s.", expectedName);
  assert(manifest.sha256(row) == hashes(inputIndex), ...
    "Analysis manifest SHA-256 mismatch for %s.", expectedName);
end
tensileRows = manifest.category == "input_tensile";
assert(nnz(tensileRows) == height(tensileInputs), ...
  "Analysis manifest must contain both registered tensile inputs.");
for inputIndex = 1:height(tensileInputs)
  expectedName = tensileInputs.name(inputIndex);
  row = tensileRows & manifest.name == expectedName;
  assert(nnz(row) == 1, ...
    "Analysis manifest tensile key is missing or duplicated: %s", ...
    expectedName);
  assert(manifest.path(row) == tensileInputs.relative_path(inputIndex), ...
    "Analysis manifest tensile path mismatch for %s.", expectedName);
  assert(manifest.sha256(row) == tensileHashes(inputIndex), ...
    "Analysis manifest tensile SHA-256 mismatch for %s.", expectedName);
end
end

function write_manuscript_readme(outputRoot, contract, catalog)
readmePath = fullfile(outputRoot, "README.md");
fileId = fopen(char(readmePath), "w", "n", "UTF-8");
assert(fileId >= 0, "Could not write README: %s", readmePath);
cleanupFile = onCleanup(@() fclose(fileId));
lines = [ ...
  "# Gr4B23271 冷变形 EBSD 综合分析"
  ""
  "本目录由统一的 MTEX 流程从六组原始 CTF 及其配对去噪数据生成。图中水平方向为棒材轴向 AD（试样 x 方向），竖直方向为面内横向/径向 TD/RD（试样 y 方向），扫描法向为 ND（试样 z 方向）。"
  ""
  "## 数据作用与分析参数"
  ""
  "原始 CTF 是定量结论的主要来源；去噪 CTF 仅用于逐状态敏感性比较，不作为独立重复。主要晶粒重构采用 2° 检测下限、15° 高角度晶界划分阈值及 minPixel=5，不进行晶界平滑。晶界统计按物理线段长度加权，同时保留 0.5°、1°、2° 和 5° 检测下限敏感性结果。"
  ""
  "KAM 采用一阶邻域和 5° 排除阈值作为主要描述，同时输出一阶/二阶邻域与 2°/5° 阈值组合的敏感性结果。GROD 轴统计仅纳入 GROD 角不低于 5° 的像素，并以反极性二阶矩张量汇总轴向；KAM、GROD 和 GOS 均为取向梯度或晶内取向离散代理量。"
  ""
  "织构采用 De la Vallee Poussin 核，半宽 5°，统一 5° 评估网格。像素加权 ODF 为主要状态描述，面积加权晶粒平均 ODF 用于权重敏感性。clipped_grid_entropy 在将负 ODF 网格值截断为零并归一至单位均值后计算；原始—去噪 ODF 距离使用单位均值谐波 ODF 的 L2 距离。"
  ""
  "c 轴作为无方向轴处理：先令 AD 分量为正，若 AD 分量在 1e-12 容差内为零，则依次令 TD/RD、ND 分量为正。由此得到的 c 轴—AD 夹角为 0–90° 锐角，不能将峰分裂、方位重分配或织构展宽简化为单一净旋转。"
  ""
  "轴向取向倾向以 AD 单轴拉伸、各滑移/孪生族最大绝对 Schmid 因子为几何假设；涵盖 basal <a>、prismatic <a>、pyramidal <a>、pyramidal <c+a> 及已验证的拉伸/压缩孪生族。未给定 CRSS，因而不据此判定实际激活顺序；Taylor 因子未计算，输出中的该列保持 NaN。"
  ""
  "## 拉伸数据纳入规则"
  ""
  "拉伸整合保留每条曲线原始 yield_status，并依据用户排除表确定是否纳入状态汇总。Rp0.2 仅汇总 yield_status=ok、数值有限且被用户接受的曲线。6.02 mm 状态的用户排除表仅接受 Y-1；Y-2 的 Rp0.2 标记为 not_reliable_no_elastic_segment_or_offset_crossing，虽其他拉伸字段仍为有限值，本流程仍依据用户排除表不将 Y-2 纳入状态汇总。因此该状态各汇总量的 n_valid=1，标准差为 0 仅表示单值汇总，不能解释为无离散性。"
  "论文面向的拉伸汇总与 EBSD 相关性仅纳入 Rp0.2_MPa、UTS_engineering_MPa 和 strain_at_UTS_uniform_elongation_percent。max_true_stress_MPa 与 true_strain_at_max_true_stress_percent 由颈缩后均匀变形公式继续换算，不能作为局部真实应力—应变，故不纳入论文面向的汇总与相关性分析。"
  ""
  "## EBSD—拉伸合并表的行粒度"
  ""
  "`08_tensile_integration/ebsd_tensile_merged.csv` 是 EBSD 指标 × 拉伸曲线 × 拉伸量的笛卡尔长表，其自然键为 `(sample, ebsd_metric, tensile_sample, tensile_repeat, tensile_metric)`。同一 EBSD 值会随不同拉伸曲线和拉伸量重复，同一 tensile_value 以及 tensile_mean、tensile_sd、tensile_n_valid 也会随不同 EBSD 指标重复；这些重复值用于保留每个指标组合的可追溯性，不能作为独立实验观测或额外重复数。"
  ""
  "`included_in_aggregate` 仅表示该条拉伸曲线的该项数值是否进入相应条件汇总；条件级均值、标准差和有效重复数以 `tensile_condition_aggregates.csv` 为准。相关性始终使用六个条件级状态，不按 `ebsd_tensile_merged.csv` 的行数扩增样本量。"
  ""
  "## 结果的证据属性"
  ""
  "| 属性 | 本流程中的结果 | 解释边界 |"
  "| --- | --- | --- |"
  "| 直接测量 | CTF 中的相、取向与质量字段，以及拉伸曲线得到的应力和应变指标 | 单个扫描区域或有限拉伸重复所代表的实验结果 |"
  "| 取向派生量 | 晶粒形貌、晶界角度、织构和 c 轴分布 | 由注册重构与统计参数计算，不等同于变形机制的直接观测 |"
  "| 代理指标 | KAM、GROD、GOS 及低角度晶界密度 | 反映取向梯度或亚结构特征，不是总位错密度的直接测量 |"
  "| 推断 | 位错储存、动态回复、滑移或孪生对加工硬化的可能贡献 | 必须同时结合力学趋势、显微组织及文献；EBSD 与 Schmid/Taylor 结果不能唯一识别实际激活机制 |"
  ""
  "## 原始—去噪比较与拉伸相关性"
  ""
  "`07_raw_denoised_comparison` 给出每个注册汇总量的配对差值、相对差值、秩次一致性和相邻变形状态趋势一致性。若去噪改变相邻趋势方向或使带符号量跨越零值，结果会被标记为可能改变科学解释；原始数值仍保持主要地位。"
  ""
  "`08_tensile_integration` 中的相关性采用六个直径状态的汇总值计算描述性 Spearman 系数，并报告逐一剔除状态后的范围。像素数未被视为实验重复，不进行高参数拟合；这些相关性用于提出假设，不能作为因果关系的证明。"
  ""
  "## 拉伸相关性指标白名单"
  ""
  "原始—去噪敏感性表保留全部数值型汇总量；拉伸相关性只纳入下列预先注册的科学指标，排除 grain_count、valid_pixel_count、total_area 和总长度等扫描支持量，并避免将同一 GROD/GOS 在不同 KAM 参数行中重复计算。"
  ""
  "- 晶粒形貌：面积加权中值等效直径、面积加权中值长宽比、数目中值形状因子、面积加权长轴—AD 中值夹角。"
  "- 晶界：晶界线密度、2–15° 低角度晶界长度分数、候选孪晶界长度分数。"
  "- 晶内取向：主要参数（一阶邻域、5° 排除阈值）下的平均 KAM、平均 GROD 和面积加权平均 GOS。"
  "- 织构：像素加权 texture index、M-index、c 轴—AD 平均夹角、c 轴位于 AD 30° 内比例及绕 AD 方位合成量。"
  "- 轴向取向倾向：各滑移/孪生族面积加权平均最大绝对 Schmid 因子。"
  ""
  "## 完整输出清单"
  ""];
lines = [lines; artifact_inventory_lines(contract, catalog); ""; ...
  "以上清单包含全部注册产物及生成器写出的附加诊断表。所有注册输入的 SHA-256、MATLAB/MTEX 版本和主要参数记录于 `analysis_manifest.csv`。"];
for lineIndex = 1:numel(lines)
  fprintf(fileId, "%s\n", lines(lineIndex));
end
end

function lines = artifact_inventory_lines(contract, catalog)
lines = ["- `README.md`"; "- `analysis_manifest.csv`"];
artifactFields = string(fieldnames(contract.artifacts));
for directoryIndex = 1:numel(contract.directories)
  directory = contract.directories(directoryIndex);
  lines(end+1,1) = ""; %#ok<AGROW>
  lines(end+1,1) = "### `" + directory + "`"; %#ok<AGROW>
  if directory == "01_standard_maps"
    for catalogIndex = 1:height(catalog)
      artifact = catalog.sample(catalogIndex) + "_" + ...
        catalog.variant(catalogIndex) + "_maps.png";
      lines(end+1,1) = "- `" + directory + "/" + artifact + "`"; %#ok<AGROW>
    end
  else
    fieldName = artifactFields(directoryIndex + 1);
    for artifact = string(contract.artifacts.(fieldName))'
      lines(end+1,1) = "- `" + directory + "/" + ...
        artifact + "`"; %#ok<AGROW>
    end
  end
  extras = additional_artifacts(directory);
  for artifact = extras'
    lines(end+1,1) = "- `" + directory + "/" + artifact + ...
      "`（附加诊断）"; %#ok<AGROW>
  end
end
end

function artifacts = additional_artifacts(directory)
switch directory
  case "02_grain_morphology"
    artifacts = "grain_morphology_quantiles.csv";
  case "03_boundaries"
    artifacts = "boundary_detection_sensitivity.csv";
  otherwise
    artifacts = strings(0,1);
end
end

function verify_output_contract(outputRoot, contract, catalogHeight)
for directory = contract.directories'
  assert(isfolder(fullfile(outputRoot, directory)), ...
    "Missing output directory: %s", directory);
end
for rootArtifact = contract.artifacts.root'
  assert(isfile(fullfile(outputRoot, rootArtifact)), ...
    "Missing root artifact: %s", rootArtifact);
end
artifactFields = string(fieldnames(contract.artifacts));
for fieldIndex = 2:numel(artifactFields)
  fieldName = artifactFields(fieldIndex);
  directoryIndex = fieldIndex - 1;
  directory = contract.directories(directoryIndex);
  artifacts = contract.artifacts.(fieldName);
  if directory == "01_standard_maps"
    maps = dir(fullfile(outputRoot, directory, "*_maps.png"));
    assert(numel(maps) == catalogHeight, ...
      "Expected one standard-map panel per scan variant.");
    continue
  end
  for artifact = string(artifacts)'
    assert(isfile(fullfile(outputRoot, directory, artifact)), ...
      "Missing registered artifact: %s/%s", directory, artifact);
  end
end
end
