function state = initChemistryV2(state)
% initChemistryV2 Initialize chemistry runtime data for mass-fraction transport.

state.chem = struct();

setupPath = getChemistrySetupPath(state);
cfg = parseChemistrySetup(setupPath);
cfg.database_path = resolvePhreeqcDatabasePathV2(cfg.database_path);

aqNames = string(cfg.aqueous_names(:));
aqSuffix = string(cfg.aqueous_suffix(:));
for i = 1:numel(aqSuffix)
    aqSuffix(i) = sanitize_aq_suffix(aqNames(i), aqSuffix(i));
end
nComp = numel(aqNames);
nCells = state.n;

aqFormula = strings(nComp, 1);
for i = 1:nComp
    aqFormula(i) = resolve_aq_formula(aqNames(i), aqSuffix(i));
end

molarMass = zeros(nComp, 1);
for i = 1:nComp
    molarMass(i) = phreeqcGfwKgV2(aqFormula(i), cfg.database_path);
end

[aqMassfrac, aqMolal] = resolve_aqueous_composition(cfg, molarMass);
plotLabels = aqNames;
for i = 1:nComp
    sx = strtrim(string(aqSuffix(i)));
    if strlength(sx) > 0
        plotLabels(i) = regexprep(sx, '^\s*as\s+', '', 'ignorecase');
    end
end

gasNames = string(cfg.gas_names(:));
[gasCompIdx, gasMask] = map_gases_to_components(gasNames, aqNames);

kinNames = string(cfg.kinetics_names(:));
nMin = numel(kinNames);
kinTau = cfg.mineral_tau_s(:);
if numel(kinTau) < nMin
    kinTau(end+1:nMin, 1) = 3600;
elseif numel(kinTau) > nMin
    kinTau = kinTau(1:nMin);
end

kinDensity = cfg.mineral_density_kg_m3(:);
if numel(kinDensity) < nMin
    kinDensity(end+1:nMin, 1) = 2500;
elseif numel(kinDensity) > nMin
    kinDensity = kinDensity(1:nMin);
end

kinFormula = string(cfg.mineral_formula(:));
if numel(kinFormula) < nMin
    kinFormula(end+1:nMin, 1) = "";
elseif numel(kinFormula) > nMin
    kinFormula = kinFormula(1:nMin);
end
kinMM = zeros(nMin, 1);
for i = 1:nMin
    f = strtrim(kinFormula(i));
    if strlength(f) == 0
        f = kinNames(i);
    end
    kinFormula(i) = f;
    kinMM(i) = phreeqcGfwKgV2(f, cfg.database_path);
end

species = repmat(struct( ...
    'name', "", ...
    'ppm', NaN, ...
    'type', "aq", ...
    'molal', 0, ...
    'formula', "", ...
    'tau', 0, ...
    'density', 0, ...
    'A', 0, ...
    'E', 0, ...
    'n', 1), 1, nComp);
for i = 1:nComp
    species(i).name = aqNames(i);
    species(i).molal = aqMolal(i);
    species(i).formula = aqFormula(i);
    if gasMask(i)
        species(i).type = "g";
    else
        species(i).type = "aq";
    end
end

state.mChem = nComp;
state.chemNames = string(aqNames(:)).';
state.molar_mass = molarMass(:).';
state.C = repmat(aqMassfrac(:), 1, nCells);
state.molal_init = aqMolal(:).';
state.massfrac_init = aqMassfrac(:).';
state.inlet_massfrac = state.massfrac_init;

state.chem.mode = "v2";
state.chem.enabled = true;
state.chem.initialized = false;
state.chem.setup = cfg;
state.chem.db = char(cfg.database_path);
state.chem.useIPhreeqc = true;
state.chem.partitionMode = 2;
if isfield(cfg, 'gas_partition') && ~isempty(cfg.gas_partition) && isfinite(cfg.gas_partition)
    state.chem.partitionMode = round(cfg.gas_partition);
end
switch state.chem.partitionMode
    case 1
        state.chem.useAnalyticalGasPartition = false;
        state.chem.gasPartitionModel = "phreeqc_gas_phase";
    case 2
        state.chem.useAnalyticalGasPartition = true;
        state.chem.gasPartitionModel = "legacy_partitionCoefficients";
    otherwise
        error('initChemistryV2:InvalidPartitionMode', ...
            'Unsupported gas partition mode %g. Use 1 (PHREEQC) or 2 (MATLAB).', ...
            state.chem.partitionMode);
end
state.chem.rm = [];
state.chem.iph = [];
state.chem.lastPhreeqcScript = "";
state.chem.nCells = nCells;
state.chem.prefix = 'chemistry_v2';
semiImplicitMaxIterDefault = max(1, round(get_state_numeric(state, 'chem_semi_implicit_max_iter', 3)));
semiImplicitRelaxationDefault = min(max(get_state_numeric(state, 'chem_semi_implicit_relaxation', 0.35), 0), 1);
state.chem.semiImplicitEnabled = true;
state.chem.semiImplicitMaxIter = semiImplicitMaxIterDefault;
state.chem.semiImplicitRelaxation = semiImplicitRelaxationDefault;
state.chem.semiImplicitKTolerance = 1e-3;
state.chem.semiImplicitMolesTolerance = 1e-6;
state.chem.lastSemiImplicitIterations = 0;
state.chem.lastSemiImplicitResidualK = NaN;
state.chem.lastSemiImplicitResidualMoles = NaN;
state.chem.lastSemiImplicitConverged = false;

state.chem.aqueousNames = aqNames;
state.chem.aqueousSuffix = aqSuffix;
state.chem.aqueousFormula = aqFormula;
state.chem.plotLabels = plotLabels;
state.chem.initialMolal = aqMolal(:);
state.chem.inletMolal = aqMolal(:);
state.chem.initialMassfrac = aqMassfrac(:);
state.chem.componentMolarMass = molarMass(:);

state.chem.gasNames = gasNames;
state.chem.gasComponentIndex = gasCompIdx(:);
state.chem.gasMask = gasMask(:);

state.chem.mineralNames = kinNames;
state.chem.mineralFormula = kinFormula;
state.chem.mineralTauS = kinTau(:);
state.chem.mineralDensity = kinDensity(:);
state.chem.mineralMolarMass = kinMM(:);
state.chem.kinMoles = zeros(nMin, nCells);
state.chem.scaleMineralNames = kinNames;
state.chem.scaleMolarMass = kinMM(:);
state.chem.scaleDensity = kinDensity(:);
state.chem.mineralArealMoles = zeros(nMin, nCells);
state.chem.mineralThickness = zeros(nMin, nCells);
state.chem.totalScaleThickness = zeros(nCells, 1);
state.chem.saturationIndices = zeros(nMin, nCells);
state.chem.partitionCoefficients = zeros(nComp, nCells);
state.chem.pH = 7 * ones(1, nCells);
state.delta = zeros(nCells, 1);

state.chem.N_liq = zeros(nComp, nCells);
state.chem.N_gas = zeros(nComp, nCells);

% Compatibility fields used by plots/HDF5 writers.
state.chem.phNames = cellstr(aqNames(:));
state.chem.phSuffixes = cellstr(aqSuffix(:));
state.chem.speciesAggregation = eye(nComp);
state.chem.species = species;
state.chem.elementSymbols = aqNames(:);
state.chem.elementSpeciesIndex = (1:nComp).';
state.chem.elementMolarMass = molarMass(:);
state.chem.elementInitialMass = zeros(nComp, nCells);
for i = 1:nComp
    state.chem.elementInitialMass(i, :) = state.C(i, :);
end

if state.calc_chem == 1
    state.chem.iph = create_iphreeqc_handle(cfg.database_path);
    state.iphreeqc = 1;
end
end

function [gasCompIdx, gasMask] = map_gases_to_components(gasNames, aqNames)
nGas = numel(gasNames);
nAq = numel(aqNames);
gasCompIdx = zeros(nGas, 1);
gasMask = false(nAq, 1);

aqNorm = strings(nAq, 1);
for i = 1:nAq
    aqNorm(i) = normalize_token(aqNames(i));
end

for g = 1:nGas
    raw = string(gasNames(g));
    base = regexprep(upper(strtrim(raw)), '\(G\)$', '');
    base = strtrim(base);
    switch base
        case "CO2"
            candidates = ["C(4)"; "CO2"];
        case "H2S"
            candidates = ["S(-2)"; "H2S"];
        otherwise
            candidates = [string(base); string(raw)];
    end

    idx = 0;
    for c = 1:numel(candidates)
        q = normalize_token(candidates(c));
        hit = find(aqNorm == q, 1);
        if ~isempty(hit)
            idx = hit;
            break;
        end
    end

    gasCompIdx(g) = idx;
    if idx > 0
        gasMask(idx) = true;
    end
end
end

function t = normalize_token(s)
t = upper(strtrim(string(s)));
t = regexprep(t, '\s+', '');
end

function suffix = sanitize_aq_suffix(name, suffix)
nm = upper(regexprep(string(name), '\s+', ''));
sx = lower(strtrim(string(suffix)));

if contains(sx, "hco3") && nm ~= "C(4)"
    suffix = "";
    return
end
if contains(sx, "so4") && nm ~= "S(6)"
    suffix = "";
    return
end
if contains(sx, "h2s") && nm ~= "S(-2)"
    suffix = "";
    return
end
if contains(sx, "sio2") && nm ~= "SI"
    suffix = "";
    return
end
suffix = string(suffix);
end

function iph = create_iphreeqc_handle(dbPath)
iph = createIPhreeqcHandleV2(dbPath);
end

function v = get_state_numeric(state, fieldName, defaultValue)
v = defaultValue;
if isstruct(state) && isfield(state, fieldName)
    raw = state.(fieldName);
    if isnumeric(raw) && isscalar(raw) && isfinite(raw)
        v = raw;
    end
end
end

function [massfrac, molal] = resolve_aqueous_composition(cfg, molarMass)
n = numel(molarMass);
molal = zeros(n, 1);
ppm = zeros(n, 1);

srcMolal = cfg.aqueous_concentration_molal(:);
srcPpm = cfg.aqueous_concentration_ppm(:);

molal(1:min(n, numel(srcMolal))) = srcMolal(1:min(n, numel(srcMolal)));
ppm(1:min(n, numel(srcPpm))) = srcPpm(1:min(n, numel(srcPpm)));

massPerKgW = ppm * 1e-6;
molMask = molal > 0;
massPerKgW(molMask) = molal(molMask) .* molarMass(molMask);
molal(~molMask) = massPerKgW(~molMask) ./ molarMass(~molMask);

massfrac = massPerKgW ./ (1 + sum(massPerKgW));
end

function formula = resolve_aq_formula(name, suffix)
nm = upper(regexprep(string(name), '\s+', ''));
sx = lower(strtrim(string(suffix)));

if nm == "C(4)" || contains(sx, "co2") || contains(sx, "hco3")
    formula = "CO2";
    return
end
if nm == "S(6)" || contains(sx, "so4")
    formula = "SO4";
    return
end
if nm == "S(-2)" || contains(sx, "h2s")
    formula = "H2S";
    return
end
if nm == "SI" || contains(sx, "sio2")
    formula = "SiO2";
    return
end

clean = regexprep(string(name), '\([^)]*\)', '');
clean = strtrim(clean);
if strlength(clean) == 0
    clean = string(name);
end
formula = clean;
end
