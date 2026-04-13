function sched = loadWHPHistory(state)
% Load WHP.csv time series for pressure, temperature, and flow rates.

sched = struct('tSec', [], 'pPa', [], 'tempC', [], 'flowTotal', [], ...
    'flowBrine', [], 'flowSteam', [], 'p_init', [], 'source', '');

file = fullfile(state.SimDir, 'WHP.csv');
if exist(file, 'file') ~= 2
    return;
end

try
    opts = detectImportOptions(file);
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
    tbl = readtable(file, opts);
catch ME
    warning('Failed to read WHP.csv: %s', ME.message);
    return;
end

tbl.Properties.VariableNames = strrep(tbl.Properties.VariableNames, ' ', '');
dateCol = pickColumn(tbl, {'Date', 'Datetime', 'Time', 'Timestamp'});
pressureCol = pickColumn(tbl, {'Pressure', 'WHP', 'WellheadPressure', 'P'});
tempCol = pickColumn(tbl, {'Temperature', 'Temp', 'T'});
totalFlowCol = pickColumn(tbl, {'TotalFlow', 'Total', 'FlowTotal'});
brineFlowCol = pickColumn(tbl, {'BrineFlow', 'WaterFlow', 'LiquidFlow'});
steamFlowCol = pickColumn(tbl, {'SteamFlow', 'VaporFlow', 'GasFlow'});

if isempty(dateCol) || isempty(pressureCol)
    warning('WHP.csv missing required columns (Date, Pressure).');
    return;
end

dates = parseDateColumn(tbl.(dateCol));
pressure = toNumericColumn(tbl.(pressureCol));
tempC = [];
flowTotal = [];
flowBrine = [];
flowSteam = [];
if ~isempty(tempCol)
    tempC = toNumericColumn(tbl.(tempCol));
end
if ~isempty(totalFlowCol)
    flowTotal = toNumericColumn(tbl.(totalFlowCol));
end
if ~isempty(brineFlowCol)
    flowBrine = toNumericColumn(tbl.(brineFlowCol));
end
if ~isempty(steamFlowCol)
    flowSteam = toNumericColumn(tbl.(steamFlowCol));
end

if isempty(dates) || isempty(pressure)
    return;
end

if ~isdatetime(dates)
    try
        dates = datetime(dates);
    catch
        return;
    end
end

mask = ~isnat(dates) & isfinite(pressure);
dates = dates(mask);
pressure = pressure(mask);
tempC = applyMask(tempC, mask);
flowTotal = applyMask(flowTotal, mask);
flowBrine = applyMask(flowBrine, mask);
flowSteam = applyMask(flowSteam, mask);
if isempty(dates)
    return;
end

baseDate = dates(1);
tSec = seconds(dates - baseDate);
[tSec, sortIdx] = sort(tSec(:));
pressure = applyIndex(pressure, sortIdx);
tempC = applyIndex(tempC, sortIdx);
flowTotal = applyIndex(flowTotal, sortIdx);
flowBrine = applyIndex(flowBrine, sortIdx);
flowSteam = applyIndex(flowSteam, sortIdx);
[tSec, uniqIdx] = unique(tSec, 'stable');
pressure = applyIndex(pressure, uniqIdx);
tempC = applyIndex(tempC, uniqIdx);
flowTotal = applyIndex(flowTotal, uniqIdx);
flowBrine = applyIndex(flowBrine, uniqIdx);
flowSteam = applyIndex(flowSteam, uniqIdx);

sched.tSec = tSec;
% WHP.csv pressure is in barg: convert to absolute bar, then to Pa.
sched.pPa = (pressure(:) + 1) * 1e5;
sched.tempC = tempC;
sched.flowTotal = flowTotal;
sched.flowBrine = flowBrine;
sched.flowSteam = flowSteam;
sched.source = file;
end

function colName = pickColumn(tbl, candidates)
colName = '';
vars = tbl.Properties.VariableNames;
varsNorm = lower(regexprep(vars, '[^a-z0-9]', ''));
for i = 1:numel(candidates)
    candNorm = lower(regexprep(candidates{i}, '[^a-z0-9]', ''));
    idx = find(strcmp(varsNorm, candNorm), 1);
    if ~isempty(idx)
        colName = vars{idx};
        return;
    end
end
end

function dt = parseDateColumn(col)
dt = [];
if isempty(col)
    return;
end

if isdatetime(col)
    dt = col;
    return;
end

if isnumeric(col)
    try
        dt = datetime(col, 'ConvertFrom', 'datenum');
        return;
    catch
        dt = [];
    end
end

if ischar(col)
    col = cellstr(col);
end

if iscell(col) || isstring(col)
    str = string(col);
    str = strtrim(str);
    try
        dt = datetime(str, 'InputFormat', 'M/d/uuuu');
    catch
        try
            dt = datetime(str);
        catch
            dt = [];
        end
    end
end
end

function out = toNumericColumn(col)
if isempty(col)
    out = [];
    return;
end
if isnumeric(col)
    out = col;
    return;
end

if isstring(col)
    col = cellstr(col);
elseif ischar(col)
    col = cellstr(col);
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

function out = applyMask(col, mask)
out = [];
if isempty(col) || numel(col) ~= numel(mask)
    return;
end
out = col(mask);
end

function out = applyIndex(col, idx)
out = [];
if isempty(col)
    return;
end
if max(idx) > numel(col)
    return;
end
out = col(idx);
end
