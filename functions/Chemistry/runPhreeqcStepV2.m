function [NliqOut, NgasOut, out] = runPhreeqcStepV2(state, NliqIn, NgasIn, dtS, scaleS, doKinetics, geom)
% runPhreeqcStepV2 Run one PHREEQC chemistry update for all cells.

iph = state.chem.iph;
[~, n] = size(NliqIn);
dtCell = reshape(dtS, 1, []);
if isscalar(dtCell)
    dtCell = repmat(dtCell, 1, n);
end
scaleCell = reshape(scaleS, 1, []);
if isscalar(scaleCell)
    scaleCell = repmat(scaleCell, 1, n);
end

aqNames = string(state.chem.aqueousNames(:));
aqSuffix = string(state.chem.aqueousSuffix(:));
gasNames = string(state.chem.gasNames(:));
gasCompIdx = state.chem.gasComponentIndex(:);
gasMask = state.chem.gasMask(:);
kinNames = string(state.chem.mineralNames(:));
kinTau = state.chem.mineralTauS(:);
kinMolesIn = state.chem.kinMoles;

nComp = numel(aqNames);
nGas = numel(gasNames);
nMin = numel(kinNames);

tplPath = string(state.chem.setup.template_path);
[~, sec] = read_ms_sections(tplPath);

molarMassKg = state.molar_mass(:);
mLiq = geom.mLiquid;
mGas = geom.mGas;
MtotIn = (NliqIn + NgasIn) .* molarMassKg;

useAnalyticalGas = state.chem.useAnalyticalGasPartition && nGas > 0;
useExplicitGas = nGas > 0 && ~useAnalyticalGas;
hasGasTemplate = isfield(sec, 'GAS_PHASE_TEMPLATE') && ~isempty(sec.GAS_PHASE_TEMPLATE);
hasGasVolume = false(1, n);
if useExplicitGas
    hasGasVolume = isfinite(geom.Vg_L) & (geom.Vg_L > 0);
end
includeGasPhase = useExplicitGas && hasGasTemplate && any(hasGasVolume);
selectedOutputIncludesGas = includeGasPhase && nGas > 0 && selected_output_requests_gas(sec.SELECTED_OUTPUT);

Kanalytic = zeros(nComp, n);
NliqPhreeqc = NliqIn;
if useAnalyticalGas
    Kanalytic = partitionCoefficients(state, geom.T_K, geom.rho_l, geom.rho_g, ...
        geom.alpha_g, geom.alpha_l);
    Kanalytic(~gasMask, :) = 0;
    [MliqPhreeqc, ~] = splitComponentMassByPartitionV2(MtotIn, Kanalytic, mLiq, mGas, gasMask);
    NliqPhreeqc = MliqPhreeqc ./ molarMassKg;
end

[~, ~, waterMassIn] = phaseMolesToLiquidMassfracV2(NliqPhreeqc, mLiq, molarMassKg);

NgasWork = NgasIn;
if useExplicitGas && gas_template_uses_fixed_volume(sec.GAS_PHASE_TEMPLATE)
    if includeGasPhase
        NgasWork = redistribute_mapped_gases_by_volume(NgasWork, geom.Vg, gasCompIdx);
    else
        NgasWork = zeros(size(NgasWork));
    end
end

baseEnv = [
    "AQ_NAMES", join_tokens(aqNames)
    "GAS_NAMES", join_tokens(gasNames)
    "KIN_NAMES", join_tokens(kinNames)
    "SI_NAMES", join_tokens(kinNames)
    "DATABASE_PATH", string(state.chem.db)
];

molal = NliqPhreeqc ./ waterMassIn;

lines = strings(0, 1);
lines(end + 1, 1) = "DATABASE " + string(state.chem.db);
lines(end + 1, 1) = "";
lines = [lines; render_plain_block(sec.STATIC, baseEnv); ""]; %#ok<AGROW>

if nMin > 0
    lines = [lines; render_rates_block(sec.RATES, baseEnv, kinNames); ""]; %#ok<AGROW>
end

for c = 1:n
    aqRows = strings(nComp, 3);
    for i = 1:nComp
        aqRows(i, 1) = aqNames(i);
        aqRows(i, 2) = fmt_num(molal(i, c));
        aqRows(i, 3) = aqSuffix(i);
    end

    gasVolL = geom.Vg_L(c);
    hasGasPhaseCell = includeGasPhase && hasGasVolume(c);

    gasRows = strings(0, 2);
    if hasGasPhaseCell
        gasRows = strings(nGas, 2);
        for g = 1:nGas
            idx = gasCompIdx(g);
            gasRows(g, 1) = gasNames(g);
            gasRows(g, 2) = fmt_num(NgasWork(idx, c));
        end
    end

    kinRows = strings(nMin, 4);
    for j = 1:nMin
        kinRows(j, 1) = kinNames(j);
        kinRows(j, 2) = "0";
        kinRows(j, 3) = fmt_num(kinMolesIn(j, c));
        kinRows(j, 4) = fmt_num(kinTau(j));
    end

    env = [
        baseEnv
        "CELL", string(c)
        "TEMP_C", fmt_num(geom.T_C(c))
        "PRESSURE_BAR", fmt_num(geom.P_bar(c))
        "PH", fmt_num(7.0)
        "DENSITY", fmt_num(geom.rho_l(c) / 1000)
        "WATER_MASS", fmt_num(waterMassIn(c))
        "DT_S", fmt_num(dtCell(c))
        "GAS_VOLUME_L", fmt_num(gasVolL)
        "VOLUME_L", fmt_num(gasVolL)
    ];

    lines = [lines; render_block(sec.SOLUTION_TEMPLATE, env, aqRows); ""]; %#ok<AGROW>

    if includeGasPhase
        if hasGasPhaseCell
            lines = [lines; render_gas_phase_block(sec.GAS_PHASE_TEMPLATE, env, gasRows); ""]; %#ok<AGROW>
        end
    end

    if doKinetics && nMin > 0
        lines = [lines; render_kinetics_block(sec.KINETICS_TEMPLATE, env, kinRows); ""]; %#ok<AGROW>
    end
end

lines = [lines; render_selected_output_block( ...
    sec.SELECTED_OUTPUT, baseEnv, selectedOutputIncludesGas); ""]; %#ok<AGROW>
lines(end + 1, 1) = "RUN_CELLS";
lines(end + 1, 1) = "    -cells 1-" + string(n);
lines(end + 1, 1) = "END";
lines(end + 1, 1) = "";

inputScript = char(strjoin(lines, newline));
errCode = iph.RunString(inputScript);
if errCode ~= 0
    errStr = iph.GetErrorString();
    dumpPath = save_failed_phreeqc_script(state, inputScript, errStr);
    if strlength(dumpPath) > 0
        error('runPhreeqcStepV2:PhreeqcFailed', ...
            'PHREEQC failed: %s\nScript saved to: %s', errStr, dumpPath);
    end
    error('runPhreeqcStepV2:PhreeqcFailed', 'PHREEQC failed: %s', errStr);
end

headings = get_selected_output_headings(iph, 1);
raw = iph.GetSelectedOutputValues(1);
data = to_numeric_matrix(raw);

rows = select_selected_output_rows_by_solution(data, headings, n);
sel = data(rows, :);

pH = extract_selected_output_scalar(sel, headings, ["pH"; "ph"], "pH");
waterMassOut = extract_selected_output_scalar( ...
    sel, headings, ["mass_H2O"; "massh2o"; "water"], "mass_H2O");
totalsMolal = extract_selected_output_totals(sel, headings, aqNames);
NliqOut = totalsMolal .* waterMassOut;

gasMoles = zeros(nGas, n);
gasWaterMoles = zeros(1, n);
if selectedOutputIncludesGas
    [gasMoles, gasWaterMoles] = extract_selected_output_gas_moles(sel, headings, gasNames);
end

si = zeros(nMin, n);
if nMin > 0
    si = extract_selected_output_si(sel, headings, kinNames);
end

kin = kinMolesIn;
if doKinetics && nMin > 0
    kinParcel = extract_selected_output_kinetics(sel, headings, kinNames);
    kin = kinMolesIn + (kinParcel - kinMolesIn) .* scaleCell;
end

NgasOut = zeros(nComp, n);
if useAnalyticalGas
    MliqPhreeqc = NliqPhreeqc .* molarMassKg;
    MliqOutChem = NliqOut .* molarMassKg;
    MtotOut = MtotIn + (MliqOutChem - MliqPhreeqc);
    [MliqOut, MgasOut] = splitComponentMassByPartitionV2(MtotOut, Kanalytic, mLiq, mGas, gasMask);
    NliqOut = MliqOut ./ molarMassKg;
    NgasOut = MgasOut ./ molarMassKg;
    Kout = Kanalytic;
else
    if selectedOutputIncludesGas
        for g = 1:nGas
            idx = gasCompIdx(g);
            NgasOut(idx, :) = NgasOut(idx, :) + gasMoles(g, :);
        end
        Kout = compute_partition_coefficients_from_selected_output( ...
            NliqOut, NgasOut, waterMassOut, gasWaterMoles, molarMassKg, gasMask);
    elseif useExplicitGas
        NgasOut = NgasWork;
        Kout = state.chem.partitionCoefficients;
    else
        Kout = zeros(nComp, n);
    end
end

out = struct();
out.pH = pH;
out.si = si;
out.kinMoles = kin;
out.partitionCoefficients = Kout;
out.waterMassKg = waterMassOut;
out.gasWaterMoles = gasWaterMoles;
out.reactionDtS = dtCell;
out.reactionScale = scaleCell;
out.script = inputScript;
end

function [hasMs, sec] = read_ms_sections(path)
sec = struct( ...
    'STATIC', strings(0, 1), ...
    'RATES', strings(0, 1), ...
    'SOLUTION_TEMPLATE', strings(0, 1), ...
    'GAS_PHASE_TEMPLATE', strings(0, 1), ...
    'KINETICS_TEMPLATE', strings(0, 1), ...
    'SELECTED_OUTPUT', strings(0, 1));
hasMs = false;
ln = splitlines(string(fileread(char(path))));
active = "";
for i = 1:numel(ln)
    line = string(ln(i));
    tok = regexp(char(line), '^\s*MS_SECTION\s+([A-Za-z0-9_]+)\s*$', 'tokens', 'once', 'ignorecase');
    if ~isempty(tok)
        name = upper(string(tok{1}));
        active = name;
        sec.(char(name)) = strings(0, 1);
        hasMs = true;
        continue
    end
    if ~isempty(regexp(char(line), '^\s*MS_END\s*$', 'once', 'ignorecase'))
        active = "";
        continue
    end
    if strlength(active) == 0
        continue
    end
    sec.(char(active))(end + 1, 1) = line; %#ok<AGROW>
end
end

function out = render_plain_block(block, env)
out = strings(0, 1);
for i = 1:numel(block)
    out(end + 1, 1) = apply_env(string(block(i)), env); %#ok<AGROW>
end
end

function out = render_selected_output_block(block, env, includeGasOutput)
if nargin < 3
    includeGasOutput = true;
end
out = strings(0, 1);
for i = 1:numel(block)
    line = apply_env(string(block(i)), env);
    if ~includeGasOutput && startsWith(lower(strtrim(line)), "-gas")
        continue
    end
    out(end + 1, 1) = line; %#ok<AGROW>
end
end

function out = render_block(block, env, aqRows)
out = strings(0, 1);
for i = 1:numel(block)
    line = string(block(i));
    if contains(line, '<AQ_NAME>') || contains(line, '<AQ_C>') || contains(line, '<SUFFIX>')
        for j = 1:size(aqRows, 1)
            s = strrep(line, '<AQ_NAME>', aqRows(j, 1));
            s = strrep(s, '<AQ_C>', aqRows(j, 2));
            s = strrep(s, '<SUFFIX>', aqRows(j, 3));
            out(end + 1, 1) = apply_env(s, env); %#ok<AGROW>
        end
        continue
    end
    out(end + 1, 1) = apply_env(line, env); %#ok<AGROW>
end
end

function out = render_kinetics_block(block, env, kinRows)
out = strings(0, 1);

iName = find(contains(block, '<KIN_NAME>'), 1, 'first');
if isempty(iName)
    out = render_plain_block(block, env);
    return
end

for i = 1:iName-1
    out(end + 1, 1) = apply_env(string(block(i)), env); %#ok<AGROW>
end

iEnd = iName;
while iEnd + 1 <= numel(block)
    nxt = string(block(iEnd + 1));
    hasKinOpt = contains(nxt, '<KIN_M0>') || contains(nxt, '<KIN_M>') || ...
        contains(nxt, '<KIN_tau>');
    if ~hasKinOpt
        break
    end
    iEnd = iEnd + 1;
end

nameTpl = string(block(iName));
optTpl = string(block(iName+1:iEnd));
for j = 1:size(kinRows, 1)
    s = strrep(nameTpl, '<KIN_NAME>', kinRows(j, 1));
    out(end + 1, 1) = apply_env(s, env); %#ok<AGROW>
    for k = 1:numel(optTpl)
        o = optTpl(k);
        o = strrep(o, '<KIN_M0>', kinRows(j, 2));
        o = strrep(o, '<KIN_M>', kinRows(j, 3));
        o = strrep(o, '<KIN_tau>', kinRows(j, 4));
        out(end + 1, 1) = apply_env(o, env); %#ok<AGROW>
    end
end

for i = iEnd+1:numel(block)
    out(end + 1, 1) = apply_env(string(block(i)), env); %#ok<AGROW>
end
end

function arr = to_numeric_matrix(values)
if isnumeric(values)
    arr = values;
    return
end

[nr, nc] = size(values);
arr = zeros(nr, nc);
for r = 1:nr
    for c = 1:nc
        v = values{r, c};
        if isnumeric(v)
            x = v;
        elseif isstring(v) || ischar(v)
            x = str2double(string(v));
        else
            x = v;
        end
        arr(r, c) = x;
    end
end
end

function headings = get_selected_output_headings(iph, userNumber)
raw = iph.GetSelectedOutputHeadings(userNumber);
if isempty(raw)
    error('runPhreeqcStepV2:MissingSelectedOutputHeadings', ...
        'Selected output %d is missing headings.', userNumber);
end
if iscell(raw)
    headings = string(raw(:)).';
elseif isstring(raw)
    headings = raw(:).';
elseif ischar(raw)
    headings = string(raw);
else
    headings = string(raw(:)).';
end
headings = strtrim(headings);
end

function rows = select_selected_output_rows_by_solution(data, headings, nCells)
solnCol = find_selected_output_heading(headings, ["soln"; "solution"], "soln");
solnId = round(data(:, solnCol));
rows = zeros(nCells, 1);
for c = 1:nCells
    idx = find(solnId == c, 1, 'last');
    if isempty(idx)
        error('runPhreeqcStepV2:MissingSelectedOutputRow', ...
            'Selected output is missing rows for solution %d.', c);
    end
    rows(c) = idx;
end
end

function values = extract_selected_output_scalar(data, headings, candidates, label)
colIdx = find_selected_output_heading(headings, candidates, label);
values = data(:, colIdx).';
end

function totalsMolal = extract_selected_output_totals(data, headings, aqNames)
aqNames = string(aqNames(:));
nComp = numel(aqNames);
nRows = size(data, 1);
totalsMolal = zeros(nComp, nRows);
for i = 1:nComp
    rawName = aqNames(i);
    validName = string(matlab.lang.makeValidName(char(rawName)));
    colIdx = find_selected_output_heading(headings, ...
        ["tot_" + rawName; ...
         "tot_" + validName; ...
         rawName + "(mol/kgw)"; ...
         validName + "(mol/kgw)"; ...
         rawName; ...
         validName], ...
        rawName + " total");
    totalsMolal(i, :) = data(:, colIdx).';
end
end

function [gasMoles, gasWaterMoles] = extract_selected_output_gas_moles(data, headings, gasNames)
gasNames = string(gasNames(:));
nGas = numel(gasNames);
nRows = size(data, 1);
gasMoles = zeros(nGas, nRows);
for g = 1:nGas
    rawName = gasNames(g);
    validName = string(matlab.lang.makeValidName(char(rawName)));
    colIdx = find_selected_output_heading(headings, ...
        ["g_" + rawName; "g_" + validName; rawName; validName], ...
        "g_" + rawName);
    gasMoles(g, :) = data(:, colIdx).';
end
gasWaterMoles = extract_selected_output_scalar(data, headings, ...
    ["g_H2O(g)"; "g_H2Og"; "H2O(g)"; "H2Og"], "g_H2O(g)");
end

function si = extract_selected_output_si(data, headings, kinNames)
kinNames = string(kinNames(:));
nKin = numel(kinNames);
nRows = size(data, 1);
si = zeros(nKin, nRows);
for i = 1:nKin
    rawName = kinNames(i);
    validName = string(matlab.lang.makeValidName(char(rawName)));
    colIdx = find_selected_output_heading(headings, ...
        ["si_" + rawName; "si_" + validName; rawName; validName], ...
        "si_" + rawName);
    si(i, :) = data(:, colIdx).';
end
end

function kin = extract_selected_output_kinetics(data, headings, kinNames)
kinNames = string(kinNames(:));
nKin = numel(kinNames);
nRows = size(data, 1);
kin = zeros(nKin, nRows);
for i = 1:nKin
    rawName = kinNames(i);
    validName = string(matlab.lang.makeValidName(char(rawName)));
    colIdx = find_selected_output_heading(headings, ...
        ["k_" + rawName; "k_" + validName; rawName; validName], ...
        "k_" + rawName);
    kin(i, :) = data(:, colIdx).';
end
end

function K = compute_partition_coefficients_from_selected_output( ...
    Nliq, Ngas, waterMassKg, gasWaterMoles, molarMassKg, gasMask)

molarMassKg = molarMassKg(:);
nComp = size(Nliq, 1);
nCells = size(Nliq, 2);
mwWaterKg = 0.01801528;

liqWaterMoles = waterMassKg ./ mwWaterKg;
liqTotalMoles = liqWaterMoles + sum(Nliq, 1);
liqMassfrac = zeros(nComp, nCells);
validLiq = isfinite(liqTotalMoles) & (liqTotalMoles > 0);
if any(validLiq)
    xLiq = zeros(nComp, nCells);
    xWater = zeros(1, nCells);
    xLiq(:, validLiq) = Nliq(:, validLiq) ./ liqTotalMoles(validLiq);
    xWater(validLiq) = liqWaterMoles(validLiq) ./ liqTotalMoles(validLiq);
    meanMassLiq = sum(xLiq(:, validLiq) .* molarMassKg, 1) + xWater(validLiq) .* mwWaterKg;
    liqMassfrac(:, validLiq) = xLiq(:, validLiq) .* molarMassKg ./ meanMassLiq;
end

gasTotalMoles = gasWaterMoles + sum(Ngas, 1);
gasMassfrac = zeros(nComp, nCells);
validGas = isfinite(gasTotalMoles) & (gasTotalMoles > 0);
if any(validGas)
    yGas = zeros(nComp, nCells);
    yWater = zeros(1, nCells);
    yGas(:, validGas) = Ngas(:, validGas) ./ gasTotalMoles(validGas);
    yWater(validGas) = gasWaterMoles(validGas) ./ gasTotalMoles(validGas);
    meanMassGas = sum(yGas(:, validGas) .* molarMassKg, 1) + yWater(validGas) .* mwWaterKg;
    gasMassfrac(:, validGas) = yGas(:, validGas) .* molarMassKg ./ meanMassGas;
end

K = zeros(nComp, nCells);
gasMaskMat = repmat(gasMask(:), 1, nCells);
validK = gasMaskMat & isfinite(liqMassfrac) & (liqMassfrac > 0) & isfinite(gasMassfrac) & (gasMassfrac >= 0);
K(validK) = gasMassfrac(validK) ./ liqMassfrac(validK);
end

function idx = find_selected_output_heading(headings, candidates, label)
if nargin < 3 || strlength(string(label)) == 0
    label = "required";
end
candidates = string(candidates(:));
idx = 0;
for i = 1:numel(candidates)
    match = find(strcmpi(headings, candidates(i)), 1);
    if ~isempty(match)
        idx = match;
        return
    end
end

normHeadings = normalize_selected_output_heading(headings);
for i = 1:numel(candidates)
    normCandidate = normalize_selected_output_heading(candidates(i));
    match = find(strcmp(normHeadings, normCandidate), 1);
    if ~isempty(match)
        idx = match;
        return
    end
end

error('runPhreeqcStepV2:MissingSelectedOutputColumn', ...
    'Selected output is missing %s column.', label);
end

function token = normalize_selected_output_heading(token)
token = lower(strtrim(string(token)));
token = regexprep(token, '\(mol/kgw\)$', '');
token = regexprep(token, '^(totals_|total_|tot_|g_|gas_|si_|k_|dk_)', '');
token = regexprep(token, '[^a-z0-9]+', '');
end

function out = render_rates_block(block, env, kinNames)
out = strings(0, 1);
block = string(block(:));
kinNames = string(kinNames(:));

i = 1;
while i <= numel(block)
    line = string(block(i));
    if startsWith(strtrim(lower(line)), "rates")
        j = i;
        while j <= numel(block)
            if strcmpi(strtrim(string(block(j))), "-end")
                break
            end
            j = j + 1;
        end

        if j <= numel(block)
            seg = string(block(i:j));
            if any(contains(seg, "<KIN_NAME>"))
                for k = 1:numel(kinNames)
                    env2 = [env; "KIN_NAME", string(kinNames(k))];
                    for s = 1:numel(seg)
                        out(end + 1, 1) = apply_env(seg(s), env2); %#ok<AGROW>
                    end
                    out(end + 1, 1) = ""; %#ok<AGROW>
                end
                i = j + 1;
                continue
            end
        end
    end

    out(end + 1, 1) = apply_env(line, env); %#ok<AGROW>
    i = i + 1;
end
end

function out = render_gas_phase_block(block, env, gasRows)
out = strings(0, 1);
for i = 1:numel(block)
    line = string(block(i));
    if contains(line, '<GAS_NAME>') || contains(line, '<GAS_MOLES>')
        for g = 1:size(gasRows, 1)
            s = line;
            s = strrep(s, '<GAS_NAME>', gasRows(g, 1));
            s = strrep(s, '<GAS_MOLES>', gasRows(g, 2));
            out(end + 1, 1) = apply_env(s, env); %#ok<AGROW>
        end
        continue
    end
    out(end + 1, 1) = apply_env(line, env); %#ok<AGROW>
end
end

function tf = gas_template_uses_fixed_volume(block)
tf = false;
for i = 1:numel(block)
    line = lower(strtrim(string(block(i))));
    if startsWith(line, "-fixed_volume")
        if contains(line, "false")
            tf = false;
        else
            tf = true;
        end
        return
    end
end
end

function tf = selected_output_requests_gas(block)
tf = false;
for i = 1:numel(block)
    line = lower(strtrim(string(block(i))));
    if startsWith(line, "-gas")
        tf = true;
        return
    end
end
end

function out = redistribute_mapped_gases_by_volume(NgasIn, Vg, gasCompIdx)
out = NgasIn;
[m, ~] = size(out);
w = Vg(:).';
sumW = sum(w);
if ~(isfinite(sumW) && sumW > 0)
    out(:) = 0;
    return
end

w = w / sumW;

mapped = unique(gasCompIdx(:));
mapped = mapped(mapped >= 1 & mapped <= m);
for k = 1:numel(mapped)
    idx = mapped(k);
    totalMol = sum(out(idx, :));
    out(idx, :) = totalMol .* w;
end
end

function s = apply_env(s, env)
for i = 1:size(env, 1)
    s = strrep(s, "<" + env(i, 1) + ">", env(i, 2));
end
end

function s = join_tokens(v)
v = string(v(:));
v = v(strlength(strtrim(v)) > 0);
if isempty(v)
    s = "";
else
    s = strjoin(cellstr(v.'), ' ');
end
end

function s = fmt_num(v)
s = string(sprintf('%.16g', v));
end

function write_text_file(path, txt)
fid = fopen(char(path), 'w');
if fid < 0
    error('runPhreeqcStepV2:CannotWriteFile', 'Unable to open file for writing: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, txt, 'char');
end

function dumpPath = save_failed_phreeqc_script(state, scriptText, errText)
dumpPath = "";
try
    outDir = string(state.SimDir);
    stamp = string(datestr(now, 'yyyymmdd_HHMMSS_FFF'));
    baseName = "phreeqc_failed_" + stamp;
    dumpPath = string(fullfile(char(outDir), char(baseName + ".phr")));

    fid = fopen(char(dumpPath), 'w');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', scriptText);

    if strlength(string(errText)) > 0
        errPath = fullfile(char(outDir), char(baseName + ".err.txt"));
        fidErr = fopen(errPath, 'w');
        cleanerErr = onCleanup(@() fclose(fidErr)); %#ok<NASGU>
        fprintf(fidErr, '%s\n', char(string(errText)));
    end
catch
    dumpPath = "";
end
end
