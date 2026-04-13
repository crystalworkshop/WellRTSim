function viewerData = rv_loadResultsH5(resultsFile)
% rv_loadResultsH5 Gather metadata required to browse simulation results.

    arguments
        resultsFile {mustBeTextScalar}
    end

    resultsFile = char(resultsFile);
    if ~isfile(resultsFile)
        error('rv_loadResultsH5:FileNotFound', ...
            'Results file not found: %s', resultsFile);
    end

    viewerData = struct();
    viewerData.file = resultsFile;

    try
        viewerData.x = h5read(resultsFile, '/meta/x');
    catch ME
        error('rv_loadResultsH5:MissingMeta', ...
            'Could not read /meta/x from %s (%s)', resultsFile, ME.message);
    end
    viewerData.x = viewerData.x(:);

    try
        viewerData.Dp0 = h5read(resultsFile, '/meta/Dp0');
        viewerData.Dp0 = viewerData.Dp0(:);
    catch
        viewerData.Dp0 = [];
    end

    try
        info = h5info(resultsFile, '/profiles/time_days');
        rawTimes = h5read(resultsFile, '/profiles/time_days');
    catch ME
        error('rv_loadResultsH5:MissingProfiles', ...
            'Results file is missing /profiles/time_days (%s)', ME.message);
    end

    rawTimes = rawTimes(:);
    nTimes = max(info.Dataspace.Size);
    if nTimes <= 0
        nTimes = numel(rawTimes);
    end
    if numel(rawTimes) > nTimes
        rawTimes = rawTimes(1:nTimes);
    end
    viewerData.times = rawTimes(:).';
    if isempty(viewerData.times)
        error('rv_loadResultsH5:NoSnapshots', ...
            'No profile snapshots were found in %s', resultsFile);
    end

    viewerData.groupNames = arrayfun(@(t) sprintf('/profiles/%.6f', t), ...
        viewerData.times, 'UniformOutput', false);

    sampleIdx = findProfileGroup(resultsFile, viewerData.groupNames);
    sampleGroup = viewerData.groupNames{sampleIdx};

    [siSets, elementSets, scaleFracSets] = detectChemistryDatasets(resultsFile, sampleGroup);
    viewerData.siDatasets = siSets;
    viewerData.elementDatasets = elementSets;
    viewerData.scaleFractionDatasets = scaleFracSets;
    viewerData.siNames = stripPrefix(siSets, 'SI_');
    viewerData.elementNames = stripPrefix(elementSets, 'element_');
    viewerData.scaleMineralNames = stripPrefix(scaleFracSets, 'scaleFrac_');
    if isempty(viewerData.scaleMineralNames) && ~isempty(viewerData.siNames)
        viewerData.scaleMineralNames = viewerData.siNames;
    end
    simDir = fileparts(resultsFile);
    viewerData.scaleDensity = loadScaleDensities(resultsFile, simDir, viewerData.scaleMineralNames);
    viewerData.scaleWallArea = computeWallArea(viewerData.x, viewerData.Dp0);

    [viewerData.timeUnit, viewerData.timeUnitLabel, viewerData.timeScale] = ...
        detectTimeUnit(resultsFile);

    [viewerData.pressureUnitLabel, viewerData.pressureUnitScaleToPa] = ...
        detectPressureUnit(resultsFile, sampleGroup);
    samplePressure = tryReadDataset(resultsFile, sampleGroup, 'P', numel(viewerData.x));

    viewerData.measuredPressure = loadMeasuredCsv(simDir, 'Pressure.csv', viewerData.x, ...
        {'DepthP', 'Depth', 'Depth_m', 'DepthM'}, {'PressureMPa', 'Pressure', 'P'}, 'pressure');
    viewerData.measuredPressure = convertMeasuredPressure( ...
        viewerData.measuredPressure, viewerData.pressureUnitLabel, viewerData.x, samplePressure);

    viewerData.measuredTemperature = loadMeasuredCsv(simDir, 'Temperature.csv', viewerData.x, ...
        {'DepthT', 'Depth', 'Depth_m', 'DepthM'}, {'Temperature', 'Temp', 'T'}, 'temperature');
    viewerData.measuredVelocity = loadMeasuredMultiCsv(simDir, 'Velocity.csv', viewerData.x, ...
        {'V_L', 'V_v'}, {'DepthV', 'Depth', 'Depth_m', 'DepthM'});
    viewerData.measuredDensity = loadMeasuredMultiCsv(simDir, 'Density.csv', viewerData.x, ...
        {'rho_l', 'rho_v'}, {'DepthD', 'Depth', 'Depth_m', 'DepthM'});

    nDepth = numel(viewerData.x);
    nProfileTimes = numel(viewerData.times);
    nScaleMinerals = numel(scaleFracSets);
    scaleThickness = zeros(nDepth, nProfileTimes);
    if nScaleMinerals > 0
        scaleFractions = zeros(nScaleMinerals, nDepth, nProfileTimes);
    else
        scaleFractions = [];
    end

    for tIdx = 1:nProfileTimes
        grp = viewerData.groupNames{tIdx};
        Dp = tryReadDataset(resultsFile, grp, 'Dp', nDepth);
        if isempty(Dp) || isempty(viewerData.Dp0)
            scaleThickness(:, tIdx) = 0;
        else
            Dp = alignLength(Dp(:), nDepth, 'replicate');
            DpInit = alignLength(viewerData.Dp0(:), nDepth, 'replicate');
            scaleThickness(:, tIdx) = max(0, 0.5 * (DpInit - Dp));
        end

        for mIdx = 1:nScaleMinerals
            data = tryReadDataset(resultsFile, grp, scaleFracSets{mIdx}, nDepth);
            if isempty(data)
                data = zeros(nDepth, 1);
            else
                data = alignLength(data(:), nDepth, 'zeros');
            end
            scaleFractions(mIdx, :, tIdx) = data(:).';
        end
    end

    viewerData.scaleThickness = scaleThickness;
    viewerData.scaleFractions = scaleFractions;
    viewerData.maxScaleThicknessPerDepth = max(scaleThickness, [], 2);
    if any(viewerData.maxScaleThicknessPerDepth > 0)
        [~, viewerData.defaultDepthIndex] = max(viewerData.maxScaleThicknessPerDepth);
    else
        viewerData.defaultDepthIndex = max(1, round(nDepth / 2));
    end
end

function idx = findProfileGroup(resultsFile, grpNames)
    for idx = 1:numel(grpNames)
        try
            h5info(resultsFile, grpNames{idx});
            return;
        catch
        end
    end
    error('rv_loadResultsH5:MissingGroups', ...
        'Profile groups not found in %s', resultsFile);
end

function [siNames, elementNames, scaleNames] = detectChemistryDatasets(file, groupName)
    siNames = {};
    elementNames = {};
    scaleNames = {};
    try
        grp = h5info(file, groupName);
        names = {grp.Datasets.Name};
    catch
        names = {};
    end
    for k = 1:numel(names)
        name = names{k};
        if startsWith(name, 'SI_')
            siNames{end+1} = name; %#ok<AGROW>
        elseif startsWith(name, 'element_')
            elementNames{end+1} = name; %#ok<AGROW>
        elseif startsWith(name, 'scaleFrac_')
            scaleNames{end+1} = name; %#ok<AGROW>
        end
    end
end

function names = stripPrefix(datasetNames, prefix)
    names = strings(size(datasetNames));
    for k = 1:numel(datasetNames)
        name = string(datasetNames{k});
        if startsWith(name, prefix)
            names(k) = extractAfter(name, strlength(prefix));
        else
            names(k) = name;
        end
    end
end

function [token, label, scale] = detectTimeUnit(resultsFile)
    token = 'd';
    label = 'days';
    scale = 1.0;

    preferred = readAttribute(resultsFile, '/', 'time_unit');
    if isempty(preferred)
        preferred = readAttribute(resultsFile, '/profiles/time_days', 'tunit');
    end
    if isempty(preferred)
        preferred = readAttribute(resultsFile, '/profiles/time_days', 'units');
    end
    if isempty(preferred)
        preferred = "days";
    end

    preferred = lower(strtrim(string(preferred)));
    switch preferred
        case {"d", "day", "days"}
            token = 'd';
            label = 'days';
            scale = 1.0;
        case {"h", "hr", "hour", "hours"}
            token = 'h';
            label = 'hours';
            scale = 24.0;
        case {"s", "sec", "second", "seconds"}
            token = 's';
            label = 'seconds';
            scale = 86400.0;
        otherwise
            token = char(preferred);
            label = char(preferred);
            scale = 1.0;
    end
end

function [label, scaleToPa] = detectPressureUnit(resultsFile, sampleGroup)
    label = readAttribute(resultsFile, '/', 'pressure_unit');
    if isempty(label)
        label = readAttribute(resultsFile, [sampleGroup '/P'], 'units');
    end
    if isempty(label)
        label = readAttribute(resultsFile, '/wellhead/pressure_MPa', 'units');
    end
    if isempty(label)
        label = "bar";
    end

    label = canonicalPressureUnit(label);
    scaleToPa = pressureUnitScaleToPa(label);

    scaleAttr = readAttribute(resultsFile, '/', 'pressure_unit_scale');
    if ~isempty(scaleAttr)
        scaleAttr = double(scaleAttr(1));
        if isfinite(scaleAttr) && scaleAttr > 0
            scaleToPa = scaleAttr;
        end
    end
end

function value = readAttribute(file, path, name)
    value = [];
    try
        info = h5info(file, path);
        idx = find(strcmp({info.Attributes.Name}, name), 1);
        if ~isempty(idx)
            value = h5readatt(file, path, name);
        end
    catch
        value = [];
    end
end

function measured = loadMeasuredCsv(simDir, fileName, xRef, depthNames, valueNames, quantityType)
    measured = struct('x', [], 'value', [], 'unit', "", 'column', "");
    if nargin < 6
        quantityType = "";
    end

    csvPath = fullfile(simDir, fileName);
    if exist(csvPath, 'file') ~= 2
        return;
    end

    x = [];
    y = [];
    valueName = "";

    try
        tbl = readtable(csvPath);
    catch
        tbl = [];
    end

    if ~isempty(tbl) && ~isempty(tbl.Properties.VariableNames)
        varNames = string(tbl.Properties.VariableNames);
        depthIdx = findColumnIndex(varNames, depthNames, "depth");
        valueIdx = findColumnIndex(varNames, valueNames, "");
        if ~isempty(depthIdx)
            x = columnToNumeric(tbl{:, depthIdx});
            if isempty(valueIdx) || valueIdx == depthIdx
                valueIdx = findFirstNumericColumn(tbl, depthIdx);
            end
            if ~isempty(valueIdx)
                y = columnToNumeric(tbl{:, valueIdx});
                valueName = varNames(valueIdx);
            end
        end
    end

    if isempty(x) || isempty(y)
        try
            raw = readmatrix(csvPath);
        catch
            raw = [];
        end
        if isempty(raw) || size(raw, 2) < 2
            return;
        end
        x = raw(:, 1);
        y = raw(:, 2);
    end

    n = min(numel(x), numel(y));
    x = x(1:n);
    y = y(1:n);
    mask = isfinite(x) & isfinite(y);
    x = x(mask);
    y = y(mask);
    if isempty(x)
        return;
    end

    if ~isempty(xRef)
        L = max(xRef(:));
        if isfinite(L) && L > 0
            x = L - x;
            x = max(min(x, L), 0);
        end
    end

    [xSorted, order] = sort(x(:));
    ySorted = y(order);
    measured.x = xSorted(:).';
    measured.value = ySorted(:).';
    measured.column = valueName;
    measured.unit = detectColumnUnit(valueName, quantityType);
end

function measured = loadMeasuredMultiCsv(simDir, fileName, xRef, defaultNames, depthNames)
    measured = struct('x', [], 'values', {{}}, 'names', strings(0));

    csvPath = fullfile(simDir, fileName);
    if exist(csvPath, 'file') ~= 2
        return;
    end

    x = [];
    data = [];
    seriesNames = strings(0);

    try
        tbl = readtable(csvPath);
    catch
        tbl = [];
    end

    if ~isempty(tbl) && ~isempty(tbl.Properties.VariableNames)
        varNames = string(tbl.Properties.VariableNames);
        depthIdx = findColumnIndex(varNames, depthNames, "depth");
        if ~isempty(depthIdx)
            x = columnToNumeric(tbl{:, depthIdx});
            dataIdx = findNumericColumnIndices(tbl, depthIdx);
            if ~isempty(dataIdx)
                data = zeros(numel(x), numel(dataIdx));
                seriesNames = strings(1, numel(dataIdx));
                for i = 1:numel(dataIdx)
                    data(:, i) = columnToNumeric(tbl{:, dataIdx(i)});
                    seriesNames(i) = varNames(dataIdx(i));
                end
            end
        end
    end

    if isempty(x) || isempty(data)
        try
            raw = readmatrix(csvPath);
        catch
            raw = [];
        end
        if isempty(raw) || size(raw, 2) < 2
            return;
        end
        x = raw(:, 1);
        data = raw(:, 2:end);
    end

    n = min(numel(x), size(data, 1));
    x = x(1:n);
    data = data(1:n, :);
    mask = isfinite(x) & any(isfinite(data), 2);
    x = x(mask);
    data = data(mask, :);
    if isempty(x) || isempty(data)
        return;
    end

    if ~isempty(xRef)
        L = max(xRef(:));
        if isfinite(L) && L > 0
            x = L - x;
            x = max(min(x, L), 0);
        end
    end

    [xSorted, order] = sort(x(:));
    data = data(order, :);
    validSeries = any(isfinite(data), 1);
    data = data(:, validSeries);
    if ~isempty(seriesNames)
        seriesNames = seriesNames(validSeries);
    end
    if isempty(data)
        return;
    end

    defaultNames = string(defaultNames);
    if isempty(defaultNames)
        defaultNames = compose("Series%d", 1:size(data, 2));
    elseif numel(defaultNames) < size(data, 2)
        extraNames = compose("Series%d", numel(defaultNames)+1:size(data, 2));
        defaultNames = [defaultNames(:).' extraNames];
    else
        defaultNames = defaultNames(:).';
    end

    if isempty(seriesNames) || all(isDefaultVarName(seriesNames))
        seriesNames = defaultNames(1:size(data, 2));
    elseif numel(seriesNames) < size(data, 2)
        seriesNames = [seriesNames(:).' defaultNames(numel(seriesNames)+1:size(data, 2))];
    end

    measured.x = xSorted(:).';
    measured.values = cell(1, size(data, 2));
    for i = 1:size(data, 2)
        measured.values{i} = data(:, i).';
    end
    measured.names = seriesNames(:).';
end

function measured = convertMeasuredPressure(measured, targetUnit, xRef, pRef)
    if isempty(measured.value)
        measured.unit = string(targetUnit);
        return;
    end

    measured.value = measured.value(:).';
    measured.x = measured.x(:).';
    targetUnit = string(targetUnit);

    rawValue = measured.value;
    rawErr = pressureProfileMismatch(rawValue, measured.x, xRef, pRef);

    if strlength(string(measured.unit)) == 0
        measured.unit = targetUnit;
        return;
    end

    sourceScale = pressureUnitScaleToPa(measured.unit);
    targetScale = pressureUnitScaleToPa(targetUnit);
    if ~isfinite(sourceScale) || ~isfinite(targetScale) || sourceScale <= 0 || targetScale <= 0
        measured.unit = targetUnit;
        return;
    end

    convertedValue = rawValue * sourceScale / targetScale;
    convertedErr = pressureProfileMismatch(convertedValue, measured.x, xRef, pRef);

    if convertedErr < rawErr
        measured.value = convertedValue;
    else
        measured.value = rawValue;
    end
    measured.unit = targetUnit;
end

function err = pressureProfileMismatch(measuredValue, measuredX, xRef, pRef)
    err = inf;
    if isempty(measuredValue) || isempty(measuredX) || isempty(xRef) || isempty(pRef)
        return;
    end

    xRef = xRef(:);
    pRef = pRef(:);
    measuredX = measuredX(:);
    measuredValue = measuredValue(:);
    n = min([numel(measuredX), numel(measuredValue), numel(xRef), numel(pRef)]);
    if n == 0
        return;
    end

    pInterp = interp1(xRef, pRef, measuredX, 'linear', 'extrap');
    mask = isfinite(measuredValue) & isfinite(pInterp);
    if ~any(mask)
        return;
    end
    delta = measuredValue(mask) - pInterp(mask);
    err = sqrt(mean(delta.^2));
end

function idx = findColumnIndex(varNames, preferred, fallbackToken)
    idx = [];
    varNames = string(varNames);
    preferred = string(preferred);

    for k = 1:numel(preferred)
        match = find(strcmpi(varNames, preferred(k)), 1);
        if ~isempty(match)
            idx = match;
            return;
        end
    end

    for k = 1:numel(preferred)
        token = preferred(k);
        if strlength(token) <= 2
            continue;
        end
        match = find(contains(lower(varNames), lower(token)), 1);
        if ~isempty(match)
            idx = match;
            return;
        end
    end

    if strlength(string(fallbackToken)) > 0
        match = find(contains(lower(varNames), lower(string(fallbackToken))), 1);
        if ~isempty(match)
            idx = match;
        end
    end
end

function idx = findFirstNumericColumn(tbl, excludeIdx)
    idx = [];
    for i = 1:width(tbl)
        if any(i == excludeIdx)
            continue;
        end
        col = columnToNumeric(tbl{:, i});
        if any(isfinite(col))
            idx = i;
            return;
        end
    end
end

function idx = findNumericColumnIndices(tbl, excludeIdx)
    idx = [];
    for i = 1:width(tbl)
        if any(i == excludeIdx)
            continue;
        end
        col = columnToNumeric(tbl{:, i});
        if any(isfinite(col))
            idx(end+1) = i; %#ok<AGROW>
        end
    end
end

function data = columnToNumeric(col)
    if isempty(col)
        data = [];
        return;
    end
    if iscell(col)
        data = str2double(string(col));
    elseif isstring(col) || ischar(col)
        data = str2double(string(col));
    else
        try
            data = double(col);
        catch
            data = str2double(string(col));
        end
    end
    data = data(:);
end

function tf = isDefaultVarName(names)
    names = string(names);
    tf = false(size(names));
    for k = 1:numel(names)
        tf(k) = ~isempty(regexp(char(names(k)), '^Var\d+$', 'once'));
    end
end

function data = tryReadDataset(file, groupName, datasetName, targetLength)
    if nargin < 4
        targetLength = [];
    end
    data = [];
    try
        data = h5read(file, sprintf('%s/%s', groupName, datasetName));
        data = data(:);
    catch
        data = [];
    end
    if ~isempty(data) && ~isempty(targetLength) && numel(data) ~= targetLength
        data = alignLength(data, targetLength, 'zeros');
    end
end

function out = alignLength(data, targetLength, mode)
    if nargin < 3
        mode = 'zeros';
    end
    data = data(:);
    n = numel(data);
    if n == targetLength
        out = data;
        return;
    end
    if n > targetLength
        out = data(1:targetLength);
        return;
    end
    if n == 0
        out = zeros(targetLength, 1);
        return;
    end
    switch lower(mode)
        case 'replicate'
            out = [data; repmat(data(end), targetLength - n, 1)];
        otherwise
            out = [data; zeros(targetLength - n, 1)];
    end
end

function unit = detectColumnUnit(varName, quantityType)
    unit = "";
    varName = lower(strtrim(string(varName)));
    quantityType = lower(strtrim(string(quantityType)));
    if strlength(varName) == 0
        return;
    end

    switch quantityType
        case "pressure"
            if contains(varName, "mpa")
                unit = "MPa";
            elseif contains(varName, "kpa")
                unit = "kPa";
            elseif contains(varName, "bar")
                unit = "bar";
            elseif contains(varName, "psi")
                unit = "psi";
            elseif contains(varName, "pa")
                unit = "Pa";
            end
    end
end

function unit = canonicalPressureUnit(unit)
    unit = strtrim(string(unit));
    switch lower(unit)
        case {"pa", "pascal", "pascals"}
            unit = "Pa";
        case {"kpa"}
            unit = "kPa";
        case {"mpa"}
            unit = "MPa";
        case {"bar", "bars"}
            unit = "bar";
        case {"psi"}
            unit = "psi";
        otherwise
            unit = string(unit);
    end
end

function scale = pressureUnitScaleToPa(unit)
    unit = canonicalPressureUnit(unit);
    switch lower(unit)
        case "pa"
            scale = 1;
        case "kpa"
            scale = 1e3;
        case "mpa"
            scale = 1e6;
        case "bar"
            scale = 1e5;
        case "psi"
            scale = 6894.757293168;
        otherwise
            scale = NaN;
    end
end

function densities = loadScaleDensities(resultsFile, simDir, mineralNames)
    nMinerals = numel(mineralNames);
    densities = 2500 * ones(nMinerals, 1);
    if nMinerals == 0
        densities = zeros(0, 1);
        return;
    end

    chemistryText = readTextDataset(resultsFile, '/inputs/chemistry_md');
    if strlength(chemistryText) == 0
        chemistryPath = fullfile(simDir, 'chemistry.md');
        if exist(chemistryPath, 'file') == 2
            chemistryText = string(fileread(chemistryPath));
        end
    end

    parsed = parseDensityList(chemistryText);
    if isempty(parsed)
        return;
    end
    if numel(parsed) >= nMinerals
        densities = parsed(1:nMinerals);
    else
        densities(1:numel(parsed)) = parsed(:);
    end
end

function wallArea = computeWallArea(x, Dp0)
    wallArea = zeros(numel(x), 1);
    if isempty(x) || isempty(Dp0)
        return;
    end

    x = x(:);
    Dp0 = alignLength(Dp0(:), numel(x), 'replicate');
    if numel(x) == 1
        dxCell = 1;
    else
        edges = zeros(numel(x) + 1, 1);
        edges(2:end-1) = 0.5 * (x(1:end-1) + x(2:end));
        edges(1) = x(1) - 0.5 * (x(2) - x(1));
        edges(end) = x(end) + 0.5 * (x(end) - x(end-1));
        dxCell = abs(diff(edges));
    end
    wallArea = pi * Dp0 .* dxCell;
end

function txt = readTextDataset(file, datasetPath)
    txt = "";
    try
        raw = h5read(file, datasetPath);
    catch
        return;
    end

    if isstring(raw)
        txt = join(raw(:), newline);
        return;
    end
    if ischar(raw)
        txt = string(raw);
        return;
    end
    if iscell(raw)
        txt = string(raw{1});
        return;
    end
    if isa(raw, 'uint8')
        txt = string(native2unicode(raw(:).', 'UTF-8'));
        return;
    end
    try
        txt = string(raw(1));
    catch
        txt = "";
    end
end

function densities = parseDensityList(chemistryText)
    densities = zeros(0, 1);
    if strlength(chemistryText) == 0
        return;
    end

    token = regexp(char(chemistryText), '(?im)^\s*density_kg_m3\s*:\s*([^\r\n#]+)', 'tokens', 'once');
    if isempty(token)
        return;
    end

    numTokens = regexp(token{1}, '[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?', 'match');
    if isempty(numTokens)
        return;
    end

    densities = zeros(numel(numTokens), 1);
    for i = 1:numel(numTokens)
        densities(i) = str2double(regexprep(numTokens{i}, '[dD]', 'e'));
    end
    densities = densities(isfinite(densities));
end
