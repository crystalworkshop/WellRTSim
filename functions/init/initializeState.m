function state = initializeState(state)
% Initialize state structure for water flow system with parameters from file
%
% Input:
%   state  - Existing state structure (may contain graphics components)
%
% Output:
%   state - Initialized state structure with all parameters

% Pull geometry information and interpolants
state = read_geom(state);

% Default number of equations
state.k = 3;

% Read parameters from file
if ~isfield(state, 'SimFile') || isempty(state.SimFile)
    state.SimFile = 'params.md';
end
paramsFile = fullfile(state.SimDir, state.SimFile);
if exist(paramsFile, 'file') ~= 2
    error('%s was not found in %s', state.SimFile, state.SimDir);
end
state.paramsFile = paramsFile;
state = readParameters(state, paramsFile);
state = normalizeInitialFlowRate(state);

% Ensure scalar defaults from parameters
state = ensureField(state, 'n', 0);
state = ensureField(state, 'maxiter', 10);
state = ensureField(state, 'epsQ', 1e-5);
state = ensureField(state, 'dt', 1.0);
state = ensureField(state, 'dt_max', inf);
state = ensureField(state, 'dt_increment', 1.01);
state = ensureField(state, 'pltf', 1);
state = ensureField(state, 'PI', 0);
state = ensureField(state, 'tunit', 's');
state = ensureField(state, 'tfin', 0);
state = ensureField(state, 't_adjust', 0);
state = ensureField(state, 't_transit', 0);
state = ensureField(state, 't_save', 0);
state = ensureField(state, 'save_csv', false);
state = ensureField(state, 'P_top', 0);
state = ensureField(state, 'H_top', 0);
state = ensureField(state, 'T_top_fin', nan);
state = ensureField(state, 'P_bot', 0);
state = ensureField(state, 'T_bot', 0);
state = ensureField(state, 'feed', 0);
state = ensureField(state, 'Q_init', 0);
state = ensureField(state, 'Qtop', state.Q_init);
state.Qtop = state.Q_init;
state = ensureField(state, 'iBC_top', 1);
state = ensureField(state, 'calc_chem', 0);
state = ensureField(state, 'eps_scale', 0);
state = ensureField(state, 'IC_switch', 1);
state = ensureField(state, 'P_unit', 'Pa');
state = ensureField(state, 'g', 9.81);
state = ensureField(state, 'tau', 3600);
state = ensureField(state, 'salt_molal', 0);
state = ensureField(state, 'db', 'phreeqc.dat');
state = ensureField(state, 'iphreeqc', 1);
state = ensureField(state, 'chemsteps', 1);
state = ensureField(state, 'chem_disp', 0);
state = ensureField(state, 'chem_source', 1);
state = ensureField(state, 'stat_chem', 0);
state = ensureField(state, 'input_units', 'ppm');
state = ensureField(state, 'results_prefix', 'results');
state = ensureField(state, 'T_surface', 20);        % surface rock temperature [°C]
state = ensureField(state, 'H_q', 0);
state = ensureField(state, 'rho_r', 0);
state = ensureField(state, 'C_r', 0);
state = ensureField(state, 'k_r', 0);
state = ensureField(state, 'rho_si', 0);
state = ensureField(state, 'rho_cast', 0);
state = ensureField(state, 'C_cast', 0);
state = ensureField(state, 'kcast', 0);

% Pressure unit handling (inputs in P_unit, simulation in Pa).
[state.pressureUnitScale, state.pressureUnitLabel] = parsePressureUnit(state.P_unit);
state.P_unit = state.pressureUnitLabel;

state.P_top = state.P_top * state.pressureUnitScale;
state.P_bot = state.P_bot * state.pressureUnitScale;

% Convert productivity index from kg/(s·P_unit) to kg/(s·Pa).
state.PI = state.PI / state.pressureUnitScale;

% Top boundary-condition mode:
%   1 - pressure from WHP.csv
%   2 - flow-rate mode
%   3 - fixed outlet pressure equal to P_top
if ~isnumeric(state.iBC_top) || isempty(state.iBC_top) || ~isfinite(state.iBC_top)
    warning('initializeState:TopBC', 'Invalid iBC_top value; using 1.');
    state.iBC_top = 1;
end
state.iBC_top = round(state.iBC_top);
if ~ismember(state.iBC_top, [1, 2, 3])
    warning('initializeState:TopBC', 'Unsupported iBC_top=%g; using 1.', state.iBC_top);
    state.iBC_top = 1;
end
state.P_top_fixed = state.P_top;

% Time settings: convert final time to seconds based on unit token.
switch lower(strtrim(char(state.tunit)))
    case {'d','day','days'}
        unitSec = 86400;
    case {'h','hr','hour','hours'}
        unitSec = 3600;
    otherwise
        unitSec = 1;
end
state.timeUnitSeconds = unitSec;
state.tfin = state.tfin * unitSec; % seconds
state.save_csv = logical(state.save_csv);
state.stat_chem = max(state.stat_chem, 0);
state.profile_save_interval = state.t_save * unitSec;
if ~isfinite(state.profile_save_interval) || state.profile_save_interval <= 0
    state.profile_save_interval = max(state.dt, 1e-9);
end
state.next_profile_save_time = 0;
state.last_profile_save_time = -inf;

% Generate the spatial grid (x measured from surface downward)
state.x = linspace(0, state.Lp, state.n)';
state.dx = state.x(2) - state.x(1);
state.nc = state.n - 1;
state.sourceTerms = zeros(state.k, 1);
state.has_feedzones = false;
state.feedzone_cells = zeros(state.n, 1, 'uint8');
state.feedzone_PI = zeros(state.n, 1);
state.feedzone_P_res = zeros(state.n, 1);
state.feedzone_H_res = zeros(state.n, 1);
state.feedzone_Qm = zeros(state.n, 1);
state.feedzones = struct('depth', [], 'depth_min', [], 'depth_max', [], 'position', [], ...
    'position_min', [], 'position_max', [], 'PI', [], 'P_res', [], 'H_res', [], 'Qm', [], ...
    'iBC', [], 'cell_idx', []);

% Wellbore interpolants: work on original x reference
state.Dp = state.wellInterpolants.innerDiameterX(state.x);
state.Dp = applyDiameterTransition(state.Dp, 2);
xb = state.Lp - state.x; % x_b = 0 at bottom
state.eps_base = state.wellInterpolants.roughnessX(xb);
state.Dp = fillMissingProfile(state.Dp, 'inner diameter');
state.eps_base = fillMissingProfile(state.eps_base, 'roughness');
state.eps = state.eps_base;
state.Dp0 = state.Dp; % retain pristine diameter profile
state = initializeGravityCache(state);

% Initialize radial rock grid and thermal matrices
state = initializeRockThermalState(state);
state = initializeHeatTransferCache(state);
% If an initial temperature profile is provided, use it to set the
% rock-wall temperature along the well (Depth=0 at top; last row=bottom).
% Load measured wellbore profiles when available (pressure/temperature)
state = loadMeasuredProfiles(state);
state = loadInitialProfiles(state);
state = loadScaleProfile(state);
state = updateWellRoughness(state);
state = loadFeedzones(state);

state.interp = true;
if state.interp
    vars = {'P_sat','Pt', 'Ht', 'T_values', 'rho_values','rho_v_values','rho_l_values', ...
        'h_v_values','h_l_values','xV_values', 'mu_values', 'alpha_values', 'beta_values','HPT'};
    load("pure_water.mat", vars{:});
    state.dens = griddedInterpolant(Ht, Pt, rho_values',"linear","nearest");
    state.dens_v = griddedInterpolant(P_sat, rho_v_values',"linear","nearest");
    state.dens_l = griddedInterpolant(P_sat, rho_l_values',"linear","nearest");
    state.h_v = griddedInterpolant(P_sat, h_v_values',"linear","nearest");
    state.h_l = griddedInterpolant(P_sat, h_l_values',"linear","nearest");
    state.Temp = griddedInterpolant(Ht, Pt, T_values',"linear","nearest");
    state.x_v = griddedInterpolant(Ht, Pt, xV_values',"linear","nearest");
    state.visc = griddedInterpolant(Ht, Pt, mu_values',"linear","nearest");
end
state.HPT = [];

fprintf('State initialization completed.\n');

state = initChemistryV2(state);

% Load wellhead schedule once during initialization
state.whpSchedule = loadWHPHistory(state);
if isempty(state.whpSchedule) || ~isstruct(state.whpSchedule)
    state.whpSchedule = struct('tSec', [], 'pPa', [], 'tempC', [], ...
        'flowTotal', [], 'flowBrine', [], 'flowSteam', [], 'p_init', []);
end
state.whpSchedule.p_init = state.P_top;
% Normalize schedule arrays to a common length
if isempty(state.whpSchedule.tSec)
    state.whpSchedule.tSec = 0;
    state.whpSchedule.pPa = state.P_top;
    state.whpSchedule.tempC = nan;
    state.whpSchedule.flowTotal = nan;
    state.whpSchedule.flowBrine = nan;
    state.whpSchedule.flowSteam = nan;
end
schedLen = numel(state.whpSchedule.tSec);
fields = {'pPa','tempC','flowTotal','flowBrine','flowSteam'};
for k = 1:numel(fields)
    fname = fields{k};
    col = state.whpSchedule.(fname);
    if isempty(col)
        state.whpSchedule.(fname) = nan(schedLen, 1);
    else
        col = col(:);
        if numel(col) < schedLen
            col(end+1:schedLen, 1) = nan;
        elseif numel(col) > schedLen
            col = col(1:schedLen);
        end
        state.whpSchedule.(fname) = col;
    end
end

% Preallocate time-series storage
state.nt_profile = max(1, floor(state.tfin / max(1e-12, state.profile_save_interval)) + 1);
state.nt_well = max(1, ceil(state.tfin / max(1e-12, state.dt)) + 1);
state.nt = state.nt_profile;

state.FlowRate = zeros(1, state.nt_well);
state.tsav = zeros(1, state.nt_well);
state.WHP = zeros(1, state.nt_well);
state.RhoTop = zeros(1, state.nt_well);
state.Quality = zeros(1, state.nt_well);
state.SteamFrac = zeros(1, state.nt_well);
state.UGasTop = zeros(1, state.nt_well);
state.ULiqTop = zeros(1, state.nt_well);
state.UMixTop = zeros(1, state.nt_well);
state.RhoGasTop = zeros(1, state.nt_well);
state.RhoLiqTop = zeros(1, state.nt_well);
state.TTop = zeros(1, state.nt_well);
state.Iterations = zeros(1, state.nt_well);
state.Tol1 = nan(1, state.nt_well);
state.ChemIterations = zeros(1, state.nt_well);
state.ChemTol1 = nan(1, state.nt_well);
state.CoutPpm = nan(state.mChem, state.nt_well);
state.stepIndex = 1;
state.tt = 0;
state.cancelFlag = false;

%% Reservoir and IO placeholders
state.Y = [];
state.Y0 = [];
state.T = [];
state.profileHistory = struct('mineralThickness', state.chem.mineralThickness);
state.P_res = zeros(1, state.n);
state.H_res = zeros(1, state.n);
state.Q_mass = zeros(1, state.n);
state.Q_v = zeros(1, state.n);
state.Q_l = zeros(1, state.n);
state.Q_v_face = zeros(1, state.n + 1);
state.Q_l_face = zeros(1, state.n + 1);
state.w_v_face = zeros(1, state.n + 1);
state.w_l_face = zeros(1, state.n + 1);
state.Q = 0;
state.Dpprev = state.Dp;
state.q_heat = zeros(1, state.n);
state.Energy = zeros(1, state.n);
state.Energy_flux = zeros(1, state.n);

state.saveCounter = 0;
state.resultsCounter = 0;
state.h5_initialized = false;
state.h5_file = '';
state.h5_well_idx = 0;
state.h5_chem_idx = 0;
state.h5_prof_idx = 0;
state.cans_dial = struct('CancelRequested', false);

return;
end

function state = readParameters(state, paramsFile)
% Read parameters from a params.md file into the state structure
%
% Inputs:
%   state - state structure
%   paramsFile - Path to the params.md file or its content
%
% Output:
%   state - Updated state structure with parameters

% Read the file
fprintf('Reading parameters from file: %s\n', paramsFile);
fileID = fopen(paramsFile, 'r');
content = fscanf(fileID, '%c', inf);
fclose(fileID);

% Parse the content line by line
lines = strsplit(content, '\n');

% Process each line
for i = 1:length(lines)
    line = strtrim(lines{i});

    % Skip empty lines, headers, and comments without parameter assignment
    if isempty(line) || startsWith(line, '#') || startsWith(line, '##') || ~contains(line, '=')
        continue;
    end
    % Split line into parameter name and value
    parts = strsplit(line, {'=', '#'});
    paramName = normalizeParameterName(parts{1});
    paramValue = strtrim(parts{2});

    % Parse value (support logicals, numbers, and tokens like 's','h','d')
    if strcmpi(paramValue, 'false')
        parsed = false;
    elseif strcmpi(paramValue, 'true')
        parsed = true;
    else
        numVal = str2double(paramValue);
        if ~isnan(numVal)
            parsed = numVal;
        else
            parsed = strtrim(paramValue); % keep as string token
        end
    end
    state.(paramName) = parsed;
end

fprintf('Parameters read successfully\n');

end

function paramName = normalizeParameterName(rawName)
paramName = strtrim(rawName);
paramName = regexprep(paramName, '[^A-Za-z0-9_]', '_');
paramName = regexprep(paramName, '_+', '_');
paramName = regexprep(paramName, '^_+|_+$', '');
if isempty(paramName)
    error('initializeState:InvalidParameterName', 'Encountered empty parameter name.');
end
if isstrprop(paramName(1), 'digit')
    paramName = ['x_' paramName];
end
end

function state = ensureField(state, fieldName, defaultValue)
if ~isfield(state, fieldName) || isempty(state.(fieldName))
    state.(fieldName) = defaultValue;
end
end

function state = normalizeInitialFlowRate(state)
% Keep backward compatibility with legacy Qtop while using Q_init internally.
hasQInit = isfield(state, 'Q_init') && ~isempty(state.Q_init);
hasQTop = isfield(state, 'Qtop') && ~isempty(state.Qtop);

if hasQInit
    if hasQTop && isnumeric(state.Q_init) && isnumeric(state.Qtop) && ...
            isfinite(state.Q_init) && isfinite(state.Qtop) && state.Q_init ~= state.Qtop
        warning('initializeState:ConflictingInitFlow', ...
            'Both Q_init and legacy Qtop are set; using Q_init=%.6g kg/s.', state.Q_init);
    end
    state.Qtop = state.Q_init;
elseif hasQTop
    state.Q_init = state.Qtop;
end
end

function D = applyDiameterTransition(D, nBlend)
% Replace step changes with a linear ramp over nBlend points.
D = D(:);
raw = D;
jumps = find(abs(diff(raw)) > 0);
blendFrac = (1:nBlend)' / (nBlend + 1);
for k = 1:numel(jumps)
    i = jumps(k);
    dL = raw(i);
    dR = raw(i + 1);
    idx = i + (1:nBlend);
    idx = idx(idx <= numel(D));
    D(idx) = dL + (dR - dL) * blendFrac(1:numel(idx));
end
end

function vec = fillMissingProfile(vec, label)
vec = vec(:);
if ~any(isnan(vec))
    return;
end

valid = find(~isnan(vec));
if isempty(valid)
    error('initializeState:InvalidGeometryProfile', ...
        'Well %s profile contains only NaN values.', label);
end

firstValid = valid(1);
if firstValid > 1
    vec(1:firstValid-1) = vec(firstValid);
end
vec = fillmissing(vec, 'previous');
end

function state = loadFeedzones(state)
% Load feedzone definitions and map them to grid cells if available.

feedFile = fullfile(state.SimDir, 'Feedzones.csv');
if exist(feedFile, 'file') ~= 2
    return;
end

[pInitX, tInitX, hasInitProfiles] = getInitPTProfiles(state);
hInitX = [];
if hasInitProfiles
    hInitX = computeEnthalpyFromPT(pInitX, tInitX, state);
end

try
    opts = detectImportOptions(feedFile);
    opts.VariableNamingRule = 'preserve';
    if isprop(opts, 'ExtraColumnsRule')
        opts.ExtraColumnsRule = 'addvars';
    end
    if isprop(opts, 'MissingRule')
        opts.MissingRule = 'fill';
    end
    if isprop(opts, 'VariableNamesLine')
        opts.VariableNamesLine = 1;
    end
    if isprop(opts, 'DataLines')
        opts.DataLines = [2 inf];
    end
    tbl = readtable(feedFile, opts);
catch ME
    warning('initializeState:FeedzoneReadFailed', '%s', ME.message);
    return;
end

tbl.Properties.VariableNames = strrep(tbl.Properties.VariableNames, ' ', '_');

depthMinCol = pickColumn(tbl, {'DepthMin', 'Depth_Min', 'DepthTop', 'TopDepth', 'DepthFrom', ...
    'DepthStart', 'Depth1'});
depthMaxCol = pickColumn(tbl, {'DepthMax', 'Depth_Max', 'DepthBottom', 'BottomDepth', 'DepthTo', ...
    'DepthEnd', 'Depth2'});
depthCol = pickColumn(tbl, {'Depth', 'MeasDepth', 'MeasuredDepth'});
piCol = pickColumn(tbl, {'PI', 'ProductivityIndex', 'ProdIndex'});
pressureCol = pickColumn(tbl, {'Pressure', 'P_res'});
enthalpyCol = pickColumn(tbl, {'Enthalpy', 'H_res'});
iBCCol = pickColumn(tbl, {'iBC', 'IBC', 'iBC_switch', 'iBCswitch', 'BC', 'FeedBC'});
massFlowCol = pickColumn(tbl, {'MassFlowRate', 'MassFlow', 'Mass_Flow', 'MassFlow_kgps', ...
    'MassFlowRate_kgps', 'MassFlowRate_kg_s'});
useFallback = isempty(piCol) || (isempty(depthCol) && (isempty(depthMinCol) || isempty(depthMaxCol))) || ...
    (~hasInitProfiles && isempty(pressureCol));
if useFallback
    [depthMin, depthMax, piVal, pressure, enthalpy, massFlow, iBC] = readFeedzonesFallback(feedFile);
else
    depthMin = [];
    depthMax = [];
    if ~isempty(depthMinCol) || ~isempty(depthMaxCol)
        if isempty(depthMinCol) && ~isempty(depthCol)
            depthMin = tableColumnToNumeric(tbl.(depthCol));
        elseif ~isempty(depthMinCol)
            depthMin = tableColumnToNumeric(tbl.(depthMinCol));
        end
        if isempty(depthMaxCol) && ~isempty(depthCol)
            depthMax = tableColumnToNumeric(tbl.(depthCol));
        elseif ~isempty(depthMaxCol)
            depthMax = tableColumnToNumeric(tbl.(depthMaxCol));
        end
    elseif ~isempty(depthCol)
        depthMin = tableColumnToNumeric(tbl.(depthCol));
        depthMax = depthMin;
    end

    if ~isempty(piCol)
        piVal = tableColumnToNumeric(tbl.(piCol));
    else
        piVal = [];
    end

    if ~isempty(pressureCol)
        pressure = tableColumnToNumeric(tbl.(pressureCol));
    else
        pressure = [];
    end

    if ~isempty(enthalpyCol)
        enthalpy = tableColumnToNumeric(tbl.(enthalpyCol));
    else
        enthalpy = [];
    end

    if ~isempty(iBCCol)
        iBC = tableColumnToNumeric(tbl.(iBCCol));
    else
        iBC = [];
    end

    massFlow = [];
    if ~isempty(massFlowCol)
        massFlow = tableColumnToNumeric(tbl.(massFlowCol));
    end
end

if isempty(depthMin) || isempty(depthMax) || isempty(piVal)
    warning('Feedzones.csv missing required columns (Depth/DepthMin/DepthMax, PI).');
    return;
end

depthMin = depthMin(:);
depthMax = depthMax(:);
piVal = piVal(:);
massFlowProvided = ~isempty(massFlow);
rowCount = min([numel(depthMin), numel(depthMax), numel(piVal)]);
depthMin = depthMin(1:rowCount);
depthMax = depthMax(1:rowCount);
piVal = piVal(1:rowCount);
pressure = resizeOptional(pressure, rowCount, NaN);
enthalpy = resizeOptional(enthalpy, rowCount, NaN);
massFlow = resizeOptional(massFlow, rowCount, 0);
iBC = resizeOptional(iBC, rowCount, NaN);

if ~isempty(iBC)
    iBC(isfinite(iBC)) = double(iBC(isfinite(iBC)) ~= 0);
end


swapMask = depthMin > depthMax;
if any(swapMask)
    tmp = depthMin(swapMask);
    depthMin(swapMask) = depthMax(swapMask);
    depthMax(swapMask) = tmp;
end

valid = isfinite(depthMin) & isfinite(depthMax) & isfinite(piVal);
if ~hasInitProfiles
    valid = valid & isfinite(pressure);
end
if ~any(valid)
    warning('Feedzones.csv has no valid feedzone rows.');
    return;
end

depthMin = depthMin(valid);
depthMax = depthMax(valid);
piVal = piVal(valid);
pressure = pressure(valid);
enthalpy = enthalpy(valid);
massFlow = massFlow(valid);
iBC = iBC(valid);
massFlow(~isfinite(massFlow)) = 0;

if ~massFlowProvided && isfield(state, 'feed') && state.feed == 3
    warning('Feedzones.csv missing MassFlowRate column; feed=3 will inject zero mass flow.');
end

depthMid = 0.5 * (depthMin + depthMax);

posMin = depthMin;
posMax = depthMax;
if isfield(state, 'x') && ~isempty(state.x)
    posMin = state.x(end) - depthMax;
    posMax = state.x(end) - depthMin;
end
if isfield(state, 'Lp') && isfinite(state.Lp)
    posMin = max(0, min(state.Lp, posMin));
    posMax = max(0, min(state.Lp, posMax));
end
posMid = 0.5 * (posMin + posMax);

pressurePaFile = [];
if ~isempty(pressure)
    pressurePaFile = normalizeFeedzonePressure(pressure, depthMid, state.pressureUnitScale);
end
enthalpyJFile = [];
if ~isempty(enthalpy)
    enthalpyJFile = enthalpy(:);
    if max(enthalpyJFile) < 1e5
        enthalpyJFile = enthalpyJFile * 1e3; % kJ/kg -> J/kg
    end
end
preferFeedzonePH = isfield(state, 'feed') && state.feed == 0;
needsPrescribed = ~preferFeedzonePH && ((isfinite(iBC) & iBC == 1) | (isnan(iBC) & ~hasInitProfiles));
if any(needsPrescribed) && ~any(isfinite(enthalpyJFile))
    warning('Feedzones.csv missing Enthalpy column required for prescribed enthalpy rows.');
    return;
end
if preferFeedzonePH && ~hasInitProfiles && any(~isfinite(enthalpyJFile))
    warning('Feedzones.csv missing Enthalpy values required for feed=0 when no initial profiles are available.');
    return;
end
pressurePa = nan(size(depthMid));
enthalpyJ = nan(size(depthMid));
piValPa = piVal / state.pressureUnitScale;

feedzone_cells = zeros(state.n, 1, 'uint8');
feedzone_PI = zeros(state.n, 1);
feedzone_P_res = zeros(state.n, 1);
feedzone_H_res = zeros(state.n, 1);
feedzone_Qm = zeros(state.n, 1);
cell_idx = nan(numel(depthMin), 1);
hasDuplicates = false;
missingPrescribedH = false;
missingInitH = false;
missingFeedzoneP = false;
missingFeedzoneH = false;
xCenters = state.x(:);
dx = state.dx;

for k = 1:numel(depthMin)
    xMin = posMin(k);
    xMax = posMax(k);
    if ~isfinite(xMin) || ~isfinite(xMax)
        continue;
    end
    [idxList, weights] = getIntervalCellWeights(xCenters, dx, xMin, xMax);
    if isempty(idxList)
        continue;
    end
    available = feedzone_cells(idxList) == 0;
    if ~any(available)
        hasDuplicates = true;
        continue;
    end
    idxUse = idxList(available);
    weightsUse = weights(available);
    if isempty(weightsUse) || ~any(isfinite(weightsUse)) || sum(weightsUse) <= 0
        weightsUse = ones(size(idxUse));
    end
    weightsUse = weightsUse(:) / sum(weightsUse);
    feedzone_cells(idxUse) = 1;
    feedzone_PI(idxUse) = piValPa(k) * weightsUse;
    feedzone_Qm(idxUse) = massFlow(k) * weightsUse;

    % Pressure: for feed=0 prefer Feedzones.csv values, otherwise keep legacy behavior.
    if preferFeedzonePH && ~isempty(pressurePaFile) && isfinite(pressurePaFile(k))
        feedzone_P_res(idxUse) = pressurePaFile(k);
        pressurePa(k) = pressurePaFile(k);
    elseif hasInitProfiles
        pCells = pInitX(idxUse);
        feedzone_P_res(idxUse) = pCells;
        wMask = isfinite(pCells) & isfinite(weightsUse);
        if any(wMask)
            pressurePa(k) = sum(pCells(wMask) .* weightsUse(wMask)) / sum(weightsUse(wMask));
        end
    elseif ~isempty(pressurePaFile) && isfinite(pressurePaFile(k))
        feedzone_P_res(idxUse) = pressurePaFile(k);
        pressurePa(k) = pressurePaFile(k);
    elseif preferFeedzonePH
        missingFeedzoneP = true;
    end

    % Enthalpy: for feed=0 prefer Feedzones.csv values, otherwise keep legacy iBC logic.
    if preferFeedzonePH
        if ~isempty(enthalpyJFile) && isfinite(enthalpyJFile(k))
            feedzone_H_res(idxUse) = enthalpyJFile(k);
            enthalpyJ(k) = enthalpyJFile(k);
        elseif hasInitProfiles
            hCells = hInitX(idxUse);
            feedzone_H_res(idxUse) = hCells;
            wMask = isfinite(hCells) & isfinite(weightsUse);
            if any(wMask)
                enthalpyJ(k) = sum(hCells(wMask) .* weightsUse(wMask)) / sum(weightsUse(wMask));
            end
            missingFeedzoneH = true;
        else
            missingFeedzoneH = true;
        end
    else
        useInitH = hasInitProfiles && (isnan(iBC(k)) || iBC(k) == 0);
        usePrescribedH = (isfinite(iBC(k)) && iBC(k) == 1) || (isnan(iBC(k)) && ~hasInitProfiles);
        if ~hasInitProfiles && isfinite(iBC(k)) && iBC(k) == 0
            useInitH = false;
            usePrescribedH = true;
            missingInitH = true;
        end
        if useInitH
            hCells = hInitX(idxUse);
            feedzone_H_res(idxUse) = hCells;
            wMask = isfinite(hCells) & isfinite(weightsUse);
            if any(wMask)
                enthalpyJ(k) = sum(hCells(wMask) .* weightsUse(wMask)) / sum(weightsUse(wMask));
            else
                useInitH = false;
                missingInitH = true;
            end
        end
        if ~useInitH && usePrescribedH
            if ~isempty(enthalpyJFile) && isfinite(enthalpyJFile(k))
                feedzone_H_res(idxUse) = enthalpyJFile(k);
                enthalpyJ(k) = enthalpyJFile(k);
            elseif hasInitProfiles
                hCells = hInitX(idxUse);
                feedzone_H_res(idxUse) = hCells;
                wMask = isfinite(hCells) & isfinite(weightsUse);
                if any(wMask)
                    enthalpyJ(k) = sum(hCells(wMask) .* weightsUse(wMask)) / sum(weightsUse(wMask));
                end
                missingPrescribedH = true;
            else
                missingPrescribedH = true;
            end
        elseif ~useInitH && ~usePrescribedH && hasInitProfiles
            hCells = hInitX(idxUse);
            feedzone_H_res(idxUse) = hCells;
            wMask = isfinite(hCells) & isfinite(weightsUse);
            if any(wMask)
                enthalpyJ(k) = sum(hCells(wMask) .* weightsUse(wMask)) / sum(weightsUse(wMask));
            end
        end
    end
    cell_idx(k) = idxUse(round(numel(idxUse) / 2));
end
if hasDuplicates
    warning('Multiple feedzones mapped to the same cell; using first occurrence per cell.');
end
if missingFeedzoneP
    warning('Feedzones.csv missing Pressure for one or more feed=0 rows; using initial profile where possible.');
end
if missingFeedzoneH
    warning('Feedzones.csv missing Enthalpy for one or more feed=0 rows; using initial profile where possible.');
end
if missingPrescribedH
    warning('Feedzones.csv missing prescribed Enthalpy for one or more iBC=1 rows; falling back where possible.');
end
if missingInitH
    warning('Initial temperature/pressure profiles missing or invalid for one or more iBC=0 rows.');
end

if hasInitProfiles
    missP = ~isfinite(pressurePa) & isfinite(posMid);
    if any(missP)
        pressurePa(missP) = interp1(xCenters, pInitX, posMid(missP), 'linear', 'extrap');
    end
    missH = ~isfinite(enthalpyJ) & isfinite(posMid);
    if any(missH)
        enthalpyJ(missH) = interp1(xCenters, hInitX, posMid(missH), 'linear', 'extrap');
    end
end

state.feedzones = struct('depth', depthMid, 'depth_min', depthMin, 'depth_max', depthMax, ...
    'position', posMid, 'position_min', posMin, 'position_max', posMax, 'PI', piValPa, ...
    'P_res', pressurePa, 'H_res', enthalpyJ, 'Qm', massFlow, 'iBC', iBC, 'cell_idx', cell_idx);
state.feedzone_cells = feedzone_cells;
state.feedzone_PI = feedzone_PI;
state.feedzone_P_res = feedzone_P_res;
state.feedzone_H_res = feedzone_H_res;
state.feedzone_Qm = feedzone_Qm;
state.has_feedzones = any(feedzone_cells ~= 0);
fprintf('Loaded %d feedzone interval(s) from %s\n', numel(depthMid), feedFile);
end

function colName = pickColumn(tbl, candidates)
colName = '';
vars = tbl.Properties.VariableNames;
varsNorm = regexprep(vars, '[^a-zA-Z0-9]', '');
for i = 1:numel(candidates)
    candNorm = regexprep(candidates{i}, '[^a-zA-Z0-9]', '');
    idx = find(strcmpi(varsNorm, candNorm), 1);
    if ~isempty(idx)
        colName = vars{idx};
        return;
    end
end
end

function [depthMin, depthMax, piVal, pressure, enthalpy, massFlow, iBC] = readFeedzonesFallback(feedFile)
% Fallback parser for Feedzones.csv when headers are missing or malformed.
depthMin = [];
depthMax = [];
piVal = [];
pressure = [];
enthalpy = [];
massFlow = [];
iBC = [];

try
    raw = readcell(feedFile);
catch
    return;
end

if isempty(raw) || size(raw, 1) < 2
    return;
end

header = raw(1, :);
data = raw(2:end, :);

idxDepthMin = matchHeaderIndex(header, {'DepthMin', 'Depth_Min', 'DepthTop', 'TopDepth', 'DepthFrom', ...
    'DepthStart', 'Depth1'});
idxDepthMax = matchHeaderIndex(header, {'DepthMax', 'Depth_Max', 'DepthBottom', 'BottomDepth', 'DepthTo', ...
    'DepthEnd', 'Depth2'});
idxDepth = matchHeaderIndex(header, {'Depth', 'MeasDepth', 'MeasuredDepth'});
idxPi = matchHeaderIndex(header, {'PI'});
idxPressure = matchHeaderIndex(header, {'Pressure', 'P_res'});
idxEnthalpy = matchHeaderIndex(header, {'Enthalpy', 'H_res'});
idxiBC = matchHeaderIndex(header, {'iBC', 'IBC', 'iBC_switch', 'iBCswitch', 'BC', 'FeedBC'});
idxMass = matchHeaderIndex(header, {'MassFlowRate', 'MassFlow', 'Mass_Flow', 'MassFlow_kgps', ...
    'MassFlowRate_kgps', 'MassFlowRate_kg_s'});

if isempty(idxDepthMin) && isempty(idxDepthMax) && isempty(idxDepth)
    if size(data, 2) >= 5
        idxDepthMin = 1;
        idxDepthMax = 2;
        idxPi = 3;
        idxPressure = 4;
        idxEnthalpy = 5;
        if size(data, 2) >= 6
            idxMass = 6;
        end
    elseif size(data, 2) >= 4
        idxDepth = 1;
        idxPi = 2;
        idxPressure = 3;
        idxEnthalpy = 4;
        if size(data, 2) >= 5
            idxMass = 5;
        end
    else
        return;
    end
end

if isempty(idxPi)
    return;
end

if ~isempty(idxDepthMin) || ~isempty(idxDepthMax)
    if isempty(idxDepthMin) && ~isempty(idxDepth)
        idxDepthMin = idxDepth;
    end
    if isempty(idxDepthMax) && ~isempty(idxDepth)
        idxDepthMax = idxDepth;
    end
    if isempty(idxDepthMin) || isempty(idxDepthMax)
        return;
    end
    depthMin = cellColumnToNumeric(data(:, idxDepthMin));
    depthMax = cellColumnToNumeric(data(:, idxDepthMax));
elseif ~isempty(idxDepth)
    depthMin = cellColumnToNumeric(data(:, idxDepth));
    depthMax = depthMin;
else
    return;
end

piVal = cellColumnToNumeric(data(:, idxPi));
if ~isempty(idxPressure)
    pressure = cellColumnToNumeric(data(:, idxPressure));
end
if ~isempty(idxEnthalpy)
    enthalpy = cellColumnToNumeric(data(:, idxEnthalpy));
end
if ~isempty(idxiBC)
    iBC = cellColumnToNumeric(data(:, idxiBC));
end
if ~isempty(idxMass)
    massFlow = cellColumnToNumeric(data(:, idxMass));
end
end

function idx = matchHeaderIndex(header, candidates)
idx = [];
if isempty(header)
    return;
end
headerStr = string(header);
headerNorm = regexprep(headerStr, '[^a-zA-Z0-9]', '');
for i = 1:numel(candidates)
    candNorm = regexprep(candidates{i}, '[^a-zA-Z0-9]', '');
    hit = find(strcmpi(headerNorm, candNorm), 1);
    if ~isempty(hit)
        idx = hit;
        return;
    end
end
end

function out = cellColumnToNumeric(col)
if isempty(col)
    out = [];
    return;
end
if isnumeric(col)
    out = col;
    return;
end
out = nan(size(col));
for i = 1:numel(col)
    val = col{i};
    if ismissing(val)
        out(i) = NaN;
    elseif isnumeric(val)
        out(i) = val;
    elseif isstring(val)
        out(i) = str2double(val);
    elseif ischar(val)
        out(i) = str2double(strtrim(val));
    else
        out(i) = NaN;
    end
end
end

function out = tableColumnToNumeric(col)
if isempty(col)
    out = [];
    return;
end
if isnumeric(col)
    out = col;
    return;
end
if iscell(col)
    out = cellColumnToNumeric(col);
    return;
end
if isstring(col)
    out = str2double(col);
    return;
end
if ischar(col)
    out = str2double(string(col));
    return;
end
try
    out = double(col);
catch
    out = nan(size(col));
end
end

function out = resizeOptional(col, n, fillVal)
if nargin < 3
    fillVal = NaN;
end
if isempty(col)
    out = repmat(fillVal, n, 1);
    return;
end
col = col(:);
if numel(col) >= n
    out = col(1:n);
else
    out = [col; repmat(fillVal, n - numel(col), 1)];
end
end

function [pInitX, tInitX, hasInit] = getInitPTProfiles(state)
% Interpolate initial pressure/temperature profiles onto the simulation grid.
pInitX = [];
tInitX = [];
hasInit = false;
if ~isfield(state, 'initProfiles') || isempty(state.initProfiles) || ~isfield(state, 'x') || isempty(state.x)
    return;
end
if ~isfield(state.initProfiles, 'pressure') || ~isfield(state.initProfiles, 'temperature')
    return;
end
initP = state.initProfiles.pressure;
initT = state.initProfiles.temperature;
if isempty(initP.position) || isempty(initP.value) || isempty(initT.position) || isempty(initT.value)
    return;
end
[pPos, pVal] = sanitizeProfile(initP.position, initP.value);
[tPos, tVal] = sanitizeProfile(initT.position, initT.value);
if isempty(pPos) || isempty(tPos)
    return;
end
x = state.x(:);
if isscalar(pPos)
    pInitX = repmat(pVal(1), size(x));
else
    pInitX = interp1(pPos, pVal, x, 'linear', 'extrap');
end
if isscalar(tPos)
    tInitX = repmat(tVal(1), size(x));
else
    tInitX = interp1(tPos, tVal, x, 'linear', 'extrap');
end
hasInit = true;
end

function [posOut, valOut] = sanitizeProfile(posIn, valIn)
% Ensure profile arrays are finite, sorted, and unique in position.
posOut = posIn(:);
valOut = valIn(:);
mask = isfinite(posOut) & isfinite(valOut);
posOut = posOut(mask);
valOut = valOut(mask);
if isempty(posOut)
    return;
end
[posOut, order] = sort(posOut);
valOut = valOut(order);
[posOut, ia] = unique(posOut, 'stable');
valOut = valOut(ia);
end

function hJ = computeEnthalpyFromPT(pPa, tC, state)
% Compute enthalpy [J/kg] from pressure [Pa] and temperature [C].
hJ = nan(size(pPa));
if isempty(pPa) || isempty(tC)
    return;
end
tK = tC + 273.15;
if isfield(state, 'HPT') && ~isempty(state.HPT)
    try
        hJ = state.HPT(pPa, tK);
        return;
    catch
        hJ = nan(size(pPa));
    end
end
try
    hJ = 1e3 * IAPWS_IF97('h_pT', pPa / 1e6, tK);
catch
    for i = 1:numel(pPa)
        try
            hJ(i) = 1e3 * IAPWS_IF97('h_pT', pPa(i) / 1e6, tK(i));
        catch
            hJ(i) = NaN;
        end
    end
end
end

function [idxList, weights] = getIntervalCellWeights(xCenters, dx, xMin, xMax)
% Determine which cell centers overlap a depth interval and their weights.
idxList = [];
weights = [];
if isempty(xCenters) || ~isfinite(dx) || dx <= 0
    return;
end
if ~isfinite(xMin) || ~isfinite(xMax)
    return;
end
if xMin > xMax
    tmp = xMin;
    xMin = xMax;
    xMax = tmp;
end
if xMax <= xMin
    [~, idx] = min(abs(xCenters - xMin));
    idxList = idx;
    weights = 1;
    return;
end
cellStart = xCenters - 0.5 * dx;
cellEnd = xCenters + 0.5 * dx;
overlap = max(0, min(cellEnd, xMax) - max(cellStart, xMin));
idxList = find(overlap > 0);
if isempty(idxList)
    [~, idx] = min(abs(xCenters - 0.5 * (xMin + xMax)));
    idxList = idx;
    weights = 1;
    return;
end
weights = overlap(idxList);
end

function [scale, label] = parsePressureUnit(unitToken)
% Map pressure unit token to Pa scale and display label.
scale = 1;
label = 'Pa';
if nargin < 1 || isempty(unitToken)
    return;
end
token = lower(strtrim(char(unitToken)));
switch token
    case {'pa', 'pascal', 'pascals'}
        scale = 1;
        label = 'Pa';
    case {'kpa'}
        scale = 1e3;
        label = 'kPa';
    case {'mpa'}
        scale = 1e6;
        label = 'MPa';
    case {'bar', 'bars'}
        scale = 1e5;
        label = 'bar';
    case {'mbar'}
        scale = 1e2;
        label = 'mbar';
    otherwise
        scale = 1;
        label = char(unitToken);
end
end

function pressurePa = normalizeFeedzonePressure(pressure, depth, unitScale)
pressure = pressure(:);
depth = depth(:);
if isempty(pressure)
    pressurePa = pressure;
    return;
end
if nargin >= 3 && isfinite(unitScale) && unitScale > 0
    pressurePa = pressure * unitScale;
    return;
end

scaleCandidates = [1, 1e3, 1e5, 1e6]; % Pa, kPa, bar, MPa
targetGrad = 1e4; % Pa/m typical hydrostatic gradient
depthMask = depth > 0 & isfinite(depth);
if any(depthMask)
    depthUse = depth(depthMask);
    pressureUse = pressure(depthMask);
    err = zeros(numel(scaleCandidates), 1);
    for k = 1:numel(scaleCandidates)
        grad = median((pressureUse * scaleCandidates(k)) ./ depthUse);
        if ~isfinite(grad) || grad <= 0
            err(k) = inf;
        else
            err(k) = abs(log10(grad) - log10(targetGrad));
        end
    end
    [~, idx] = min(err);
    scale = scaleCandidates(idx);
else
    if max(pressure) > 1e5
        scale = 1;
    else
        scale = 1e5;
    end
end
pressurePa = pressure * scale;
end

function state = initializeRockThermalState(state)
% Initialize a static rock-wall temperature profile used by the heat-loss model.

state.nc = max(0, state.nc);
xc = 0.5 * (state.x(1:end-1) + state.x(2:end));

file = fullfile(state.SimDir, 'InitTemperature.csv');
state.hasRockTemperatureData = false;
if exist(file, 'file') ~= 2
    WallTemp = state.T_surface + (state.T_bot - state.T_surface) * (1 - xc /  state.Lp)+273.15; % Kelvin
else
    try
        tbl = readtable(file);
    catch ME
        error('Failed to read %s: %s', file, ME.message);
    end
    [depth, mask] = filterDepthForSim(state, tbl.DepthT);
    if isempty(depth)
        WallTemp = state.T_surface + (state.T_bot - state.T_surface) * (1 - xc /  state.Lp)+273.15; % Kelvin
    else
        dpt = flip(xc(end) - depth);
        WallTemp = interp1(dpt, flip(tbl.Temperature(mask)), xc) + 273.15; % Kelvin
        state.hasRockTemperatureData = true;
    end
end
state.T_rock = WallTemp; % rock wall temperature at wellbore [K]
end

function state = initializeHeatTransferCache(state)
% Precompute depth-dependent heat-transfer coefficient on the simulation grid.

xNodes = state.x(:);
file = fullfile(state.SimDir, 'InitTemperature.csv');
HqNode = state.H_q * ones(size(xNodes));

if exist(file, 'file') == 2
    tbl = readtable(file);
    depthCol = pickColumn(tbl, {'DepthT', 'Depth', 'Depth_m', 'DepthM', ...
        'MeasDepth', 'MeasuredDepth'});
    hqCol = pickColumn(tbl, {'H_q', 'Hq', 'HeatTransferCoeff', ...
        'HeatTransferCoefficient'});

    if ~isempty(depthCol) && ~isempty(hqCol)
        [depth, mask] = filterDepthForSim(state, tableColumnToNumeric(tbl.(depthCol)));
        position = depthToPosition(depth, state);
        hqProfile = flip(tableColumnToNumeric(tbl.(hqCol)(mask)));
        [position, hqProfile] = sanitizeProfile(position, hqProfile);
        if isscalar(position)
            HqNode(:) = hqProfile(1);
        else
            HqNode = interp1(position, hqProfile, xNodes, 'linear', 'extrap');
        end
    end
end

state.H_q_node = HqNode(:);
state.H_q_interp = griddedInterpolant(xNodes, state.H_q_node, 'linear', 'nearest');
end

function state = initializeGravityCache(state)
% Precompute gravity inclination on the simulation grid and its faces.

xNodes = state.x(:);
[thetaNode, cthNode] = getGravityInclination(state, xNodes);

if numel(xNodes) <= 1
    xFaces = xNodes;
else
    xFaces = zeros(numel(xNodes) + 1, 1);
    xFaces(1) = xNodes(1);
    xFaces(end) = xNodes(end);
    xFaces(2:end-1) = 0.5 * (xNodes(1:end-1) + xNodes(2:end));
end
[thetaFace, cthFace] = getGravityInclination(state, xFaces);

state.gravityThetaNode = thetaNode(:);
state.gravityCthNode = cthNode(:);
state.gravityThetaFace = thetaFace(:);
state.gravityCthFace = cthFace(:);
end

function state = loadMeasuredProfiles(state)
% Read measured pressure and temperature profiles if CSV files are present.
SimDir = state.SimDir;
measured = struct();
measured.pressure = struct('position', [], 'value', [], 'units', state.pressureUnitLabel);
measured.temperature = struct('position', [], 'value', [], 'units', 'degC');
measured.velocity = struct('position', [], 'values', {{}}, 'names', strings(0,1));
measured.density = struct('position', [], 'values', {{}}, 'names', strings(0,1));

% Pressure profile (expects Depth [m] and Pressure in P_unit)
pressureFile = fullfile(SimDir, 'Pressure.csv');
if exist(pressureFile, 'file') == 2
    tbl = readtable(pressureFile);
    [depth, mask] = filterDepthForSim(state, tbl.DepthP);
    if ~isempty(depth)
        position = depthToPosition(depth, state);
        if ismember('PressureMPa', tbl.Properties.VariableNames)
            value = flip(tbl.PressureMPa(mask));
        elseif ismember('Pressure', tbl.Properties.VariableNames)
            value = flip(tbl.Pressure(mask));
        elseif ismember('P', tbl.Properties.VariableNames)
            value = flip(tbl.P(mask));
        else
            value = [];
        end
        if ~isempty(value)
            measured.pressure.position = position;
            measured.pressure.value = value;
            fprintf('Loaded measured pressure profile (%d points).\n', numel(position));
        end
    end
end

% Temperature profile (expects DepthT [m] and Temperature in °C)
tempFile = fullfile(SimDir, 'Temperature.csv');
if exist(tempFile, 'file') == 2
    tbl = readtable(tempFile);
    [depth, mask] = filterDepthForSim(state, tbl.DepthT);
    if ~isempty(depth)
        position = depthToPosition(depth, state);
        value = flip(tbl.Temperature(mask));
        measured.temperature.position = position;
        measured.temperature.value = value;
        fprintf('Loaded measured temperature profile (%d points).\n', numel(position));
    end
end

% Velocity profile (DepthV plus V_L/V_v columns)
velFile = fullfile(SimDir, 'Velocity.csv');
if exist(velFile, 'file') == 2
    tbl = readtable(velFile);
    measVel = parseDepthTable(tbl, 'DepthV', state);
    if ~isempty(measVel.position)
        measured.velocity = measVel;
        fprintf('Loaded measured velocity profile (%d points, %d series).\n', ...
            numel(measVel.position), numel(measVel.values));
    end
end

% Density profile (DepthD plus rho_l/rho_v/S_v columns)
densFile = fullfile(SimDir, 'Density.csv');
if exist(densFile, 'file') == 2
    tbl = readtable(densFile);
    measDens = parseDepthTable(tbl, 'DepthD', state);
    if ~isempty(measDens.position)
        measured.density = measDens;
        fprintf('Loaded measured density profile (%d points, %d series).\n', ...
            numel(measDens.position), numel(measDens.values));
    end
end

state.measured = measured;
end

function state = loadInitialProfiles(state)
% Read initial pressure and temperature profiles for dashed overlays.
SimDir = state.SimDir;
initProfiles = struct();
initProfiles.pressure = struct('position', [], 'value', [], 'units', 'Pa');
initProfiles.temperature = struct('position', [], 'value', [], 'units', 'degC');

% Initial pressure profile (expects Depth [m] and Pressure in P_unit unless unit hinted)
pFile = fullfile(SimDir, 'InitPressure.csv');
if exist(pFile, 'file') == 2
    try
        tbl = readtable(pFile);
    catch ME
        warning('Failed to read %s: %s', pFile, ME.message);
        tbl = [];
    end
    if ~isempty(tbl)
        depthCol = pickColumn(tbl, {'DepthP', 'DepthT', 'Depth', 'Depth_m', 'DepthM', ...
            'MeasDepth', 'MeasuredDepth'});
        pressureCol = pickColumn(tbl, {'Pressure', 'P', 'PressureMPa', 'PressurekPa', ...
            'PressureKPa', 'Pressurebar', 'PressureBar', 'PressurePa', 'PressurePA'});
        if ~isempty(depthCol) && ~isempty(pressureCol)
            [depth, mask] = filterDepthForSim(state, tbl.(depthCol));
            if ~isempty(depth)
                position = depthToPosition(depth, state);
                value = flip(tbl.(pressureCol)(mask));
                scale = pressureColumnScale(pressureCol, state.pressureUnitScale);
                initProfiles.pressure.position = position;
                initProfiles.pressure.value = value * scale;
                initProfiles.pressure.units = 'Pa';
                fprintf('Loaded init pressure profile (%d points).\n', numel(position));
            end
        end
    end
end

% Initial temperature profile (expects Depth [m] and Temperature in degC)
tFile = fullfile(SimDir, 'InitTemperature.csv');
if exist(tFile, 'file') == 2
    try
        tbl = readtable(tFile);
    catch ME
        warning('Failed to read %s: %s', tFile, ME.message);
        tbl = [];
    end
    if ~isempty(tbl)
        depthCol = pickColumn(tbl, {'DepthT', 'Depth', 'Depth_m', 'DepthM', ...
            'MeasDepth', 'MeasuredDepth'});
        tempCol = pickColumn(tbl, {'Temperature', 'Temp', 'T'});
        if ~isempty(depthCol) && ~isempty(tempCol)
            [depth, mask] = filterDepthForSim(state, tbl.(depthCol));
            if ~isempty(depth)
                position = depthToPosition(depth, state);
                value = flip(tbl.(tempCol)(mask));
                initProfiles.temperature.position = position;
                initProfiles.temperature.value = value;
                initProfiles.temperature.units = 'degC';
                fprintf('Loaded init temperature profile (%d points).\n', numel(position));
            end
        end
    end
end

state.initProfiles = initProfiles;
end

function state = loadScaleProfile(state)
% Read scale thickness profile and convert it to a scaled diameter overlay.
SimDir = state.SimDir;
scaleProfile = struct('position', [], 'diameter', [], 'units', 'm');

scaleFile = fullfile(SimDir, 'scale.csv');
if exist(scaleFile, 'file') ~= 2
    state.scaleProfile = scaleProfile;
    return;
end

try
    tbl = readtable(scaleFile);
catch ME
    warning('Failed to read %s: %s', scaleFile, ME.message);
    state.scaleProfile = scaleProfile;
    return;
end

depthCol = pickColumn(tbl, {'DepthScale', 'Depth', 'Depth_m', 'DepthM', ...
    'MeasDepth', 'MeasuredDepth'});
thicknessCol = pickColumn(tbl, {'Thickness_cm', 'ThicknessCM', 'ThicknessCm', ...
    'Thickness', 'ScaleThickness', 'ScaleThickness_cm'});
if isempty(depthCol) || isempty(thicknessCol)
    state.scaleProfile = scaleProfile;
    return;
end

[depth, mask] = filterDepthForSim(state, tableColumnToNumeric(tbl.(depthCol)));
if isempty(depth)
    state.scaleProfile = scaleProfile;
    return;
end

thicknessCm = tableColumnToNumeric(tbl.(thicknessCol));
thicknessCm = thicknessCm(mask);
position = depthToPosition(depth, state);
thicknessM = 1e-2 * flip(thicknessCm(:));
baseDiameter = interp1(state.x(:), state.Dp0(:), position(:), 'linear');
scaledDiameter = baseDiameter - 2 * thicknessM;
[position, scaledDiameter] = sanitizeProfile(position, scaledDiameter);

scaleProfile.position = position;
scaleProfile.diameter = scaledDiameter;
scaleProfile.units = 'm';
state.scaleProfile = scaleProfile;

fprintf('Loaded scaled diameter profile from scale.csv (%d points).\n', numel(position));
end

function measStruct = parseDepthTable(tbl, depthField, state)
% Convert a depth-based table into plotting-ready measured struct.
measStruct = struct('position', [], 'values', {{}}, 'names', strings(0,1));
if isempty(tbl) || ~ismember(depthField, tbl.Properties.VariableNames)
    return;
end

depth = tbl.(depthField);
if isempty(depth)
    return;
end
[depth, mask] = filterDepthForSim(state, depth);
if isempty(depth)
    return;
end
pos = depthToPosition(depth, state);

varNames = tbl.Properties.VariableNames;
varNames = setdiff(varNames, depthField, 'stable');
if isempty(varNames)
    return;
end

values = cell(1, numel(varNames));
for i = 1:numel(varNames)
    col = tbl.(varNames{i});
    if numel(col) == numel(mask)
        col = col(mask);
    end
    values{i} = flip(col(:).');
end

measStruct.position = pos(:).';
measStruct.values = values;
measStruct.names = string(varNames);
end

function [depthOut, mask] = filterDepthForSim(~, depthIn)
% Normalize depth vector for profile parsing.
depthOut = depthIn(:);
mask = true(size(depthOut));
end

function pos = depthToPosition(depthVec, state)
% Convert depth (measured from surface down) to simulation x-position (0=top).
depthVec = depthVec(:);
% if ~isfield(state, 'Lp') || ~isfinite(state.Lp) || state.Lp <= 0
%     pos = flip(depthVec(end) - depthVec);
%     return;
% end
pos = state.x(end) - depthVec;
% pos = max(0, min(state.Lp, pos));
pos = flip(pos);
end

function scale = pressureColumnScale(colName, defaultScale)
% Determine pressure scale for a column name, fall back to defaultScale.
if nargin < 2 || isempty(defaultScale) || ~isfinite(defaultScale) || defaultScale <= 0
    defaultScale = 1;
end
token = lower(strtrim(char(colName)));
if contains(token, 'mpa')
    scale = 1e6;
elseif contains(token, 'kpa')
    scale = 1e3;
elseif contains(token, 'bar')
    scale = 1e5;
elseif contains(token, 'pa')
    scale = 1;
else
    scale = defaultScale;
end
end
