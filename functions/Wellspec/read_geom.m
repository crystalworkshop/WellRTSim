%% Well Casing and Deviation Analysis
% This script reads well casing and deviation data from CSV files,
% processes them, and creates gridded interpolants for various well properties.
% Visualizations are displayed in a tabbed UI figure.

function state=read_geom(state)
SimDir = state.SimDir;
% Define the input filenames
casingFilename    = [SimDir,'CASTINGDETAIL.csv'];
deviationFilename = [SimDir,'DEVIATIONDETAIL.csv'];

% Process the data and create interpolants
state.wellInterpolants = processWellCasingData(casingFilename, deviationFilename);

% Get well data from the casing file
wellData = getWellDataFromFile(casingFilename);

% Store well data in state for later use by plotting routines
state.wellData = wellData;

% Compute total length of slotted-liner region and store in state.L_feed
maskSlot = strcmp(wellData.CasingType, 'Slotted-Liner');
segLen = wellData.BottomMeasDepth(maskSlot) - wellData.TopMeasDepth(maskSlot);
segLen = segLen(~isnan(segLen) & segLen > 0);
state.L_feed = sum(segLen);
fprintf('Total slotted-liner length (L_feed): %.2f m\n', state.L_feed);

% Determine depth range
minDepth = 0;
maxDepth = max(wellData.BottomMeasDepth);

% Determine if we have deviation data
hasDeviationData = exist(deviationFilename, 'file') == 2;

% If no deviation data is available, create vertical deviation data
if ~hasDeviationData
    fprintf('No deviation data found. Assuming vertical well.\n');
    deviationData = createVerticalDeviationData(minDepth, maxDepth);
    state.wellInterpolants.offset = griddedInterpolant(deviationData.MeasDepth, deviationData.Offset, 'linear', 'linear');
    state.wellInterpolants.angle = griddedInterpolant(deviationData.MeasDepth, deviationData.Angle*0, 'linear', 'linear');
else
    % Read deviation data from file
    deviationData = getDeviationDataFromFile(deviationFilename);
end

% Store deviation data in state for later use by plotting routines
state.deviationData = deviationData;

state.Lp=maxDepth;

fprintf('Well geometry data processed successfully.\n');
end

function wellData = getWellDataFromFile(casingFilename)
% Read the casing CSV file
opts = detectImportOptions(casingFilename);
opts.VariableNamingRule = 'preserve';
opts.Delimiter = ',';
wellData = readtable(casingFilename, opts);

% Fix variable names if they have spaces
wellData.Properties.VariableNames = strrep(wellData.Properties.VariableNames, ' ', '_');

% Normalize known casing type variants from field data imports.
if ismember('CasingType', wellData.Properties.VariableNames)
    casingTypes = string(wellData.CasingType);
    casingTypes = strtrim(erase(casingTypes, '"'));
    casingTypes(strcmpi(casingTypes, "Casting")) = "Casing";
    casingTypes(strcmpi(casingTypes, "Casing")) = "Casing";
    casingTypes(strcmpi(casingTypes, "Liner")) = "Liner";
    casingTypes(strcmpi(casingTypes, "Slotted-Liner")) = "Slotted-Liner";
    wellData.CasingType = cellstr(casingTypes);
end

% Sort by depth to ensure proper processing
wellData = sortrows(wellData, 'TopMeasDepth');
 fprintf('Reading well casting data from %s\n', casingFilename);
end

function deviationData = getDeviationDataFromFile(deviationFilename)
% Read deviation data
opts = detectImportOptions(deviationFilename);
opts.VariableNamingRule = 'preserve';
opts.Delimiter = ',';
deviationData = readtable(deviationFilename, opts);

% Fix variable names
deviationData.Properties.VariableNames = strrep(deviationData.Properties.VariableNames, ' ', '_');

% Sort by measured depth
deviationData = sortrows(deviationData, 'MeasDepth');
end

function deviationData = createVerticalDeviationData(minDepth, maxDepth)
% Create vertical deviation data (zero horizontal offset)
numPoints = 50; % Number of points for vertical profile

% Create values
measDepth = linspace(minDepth, maxDepth, numPoints)';
vertDepth = measDepth;
angle = zeros(numPoints, 1);
offset = zeros(numPoints, 1);

% Create well name and ID columns
wellName = repmat({'Unknown'}, numPoints, 1);
ID = ones(numPoints, 1);

% Create table
deviationData = table(wellName, ID, measDepth, vertDepth, angle, offset, ...
    'VariableNames', {'WellName', 'ID', 'MeasDepth', 'VertDepth', 'Angle', 'Offset'});
end

function wellInterpolants = processWellCasingData(casingFilename, deviationFilename)
% Read well data
wellData = getWellDataFromFile(casingFilename);

% Check if deviation file exists and read it if available
hasDeviationData = exist(deviationFilename, 'file') == 2;
deviationData = table();
if hasDeviationData
    fprintf('Reading well deviation data from %s\n', deviationFilename);
    deviationData = getDeviationDataFromFile(deviationFilename);
end

%% Create depth array with all transition points
minDepth = 0;
maxDepth = max(wellData.BottomMeasDepth);

% Get all unique transition depths
allDepths = unique([wellData.TopMeasDepth; wellData.BottomMeasDepth]);
allDepths = sort(allDepths);

%% Generate fine sampling for visualization (1 meter increments)
depthSampling = minDepth:1:maxDepth;

%% Add deviation data if available
if hasDeviationData
    % Create offset interpolant using measured depth as the independent variable
    offsetInterpolant = griddedInterpolant(deviationData.MeasDepth, deviationData.Offset, 'linear', 'linear');

    angleData = normalizeDeviationAngles(deviationData.Angle);

    % Create angle interpolant using measured depth as the independent variable
    angleInterpolant = griddedInterpolant(deviationData.MeasDepth, angleData, 'linear', 'linear');

    % Sample offset values at our depth sampling points
    offsetSampled = offsetInterpolant(depthSampling);
    angleSampled = angleInterpolant(depthSampling);
else
    % Use zeros if no deviation data
    offsetSampled = zeros(size(depthSampling));
    angleSampled = zeros(size(depthSampling));
    offsetInterpolant = griddedInterpolant(depthSampling, offsetSampled, 'linear', 'linear');
    angleInterpolant = griddedInterpolant(depthSampling, angleSampled, 'linear', 'linear');
end

[vertDepthSampled, gravProjSampled, gravAngleSampled] = ...
    sampleGravityGeometry(depthSampling, hasDeviationData, deviationData, angleSampled);

[innerDiameterSampled, wallThicknessSampled, perforationStatusSampled, roughnessSampled] = ...
    sampleCasingProfiles(wellData, depthSampling);

%% Create the interpolants using proper step function representation
% First, determine the breakpoints (where values change)
idBreaks = findStepBreakpoints(depthSampling, innerDiameterSampled);
thicknessBreaks = findStepBreakpoints(depthSampling, wallThicknessSampled);
perfBreaks = findStepBreakpoints(depthSampling, perforationStatusSampled);
roughnessBreaks = findStepBreakpoints(depthSampling, roughnessSampled);

% Create interpolant points with step function precision
idPoints = createStepPoints(depthSampling, innerDiameterSampled, idBreaks);
thicknessPoints = createStepPoints(depthSampling, wallThicknessSampled, thicknessBreaks);
perfPoints = createStepPoints(depthSampling, perforationStatusSampled, perfBreaks);
roughnessPoints = createStepPoints(depthSampling, roughnessSampled, roughnessBreaks);

%% Create the gridded interpolants
idInterpolant = griddedInterpolant(idPoints.depths, idPoints.values, 'nearest', 'nearest');
thicknessInterpolant = griddedInterpolant(thicknessPoints.depths, thicknessPoints.values, 'nearest', 'nearest');
perfInterpolant = griddedInterpolant(perfPoints.depths, perfPoints.values, 'nearest', 'nearest');
roughnessInterpolant = griddedInterpolant(roughnessPoints.depths, roughnessPoints.values, 'nearest', 'nearest');

%% Offset and angle interpolants already created above

%% Save interpolants to MAT file
wellInterpolants = struct();
wellInterpolants.innerDiameter = idInterpolant;
wellInterpolants.wallThickness = thicknessInterpolant;
% Store top-referenced perforation status as a separate field
wellInterpolants.perforationStatusTop = perfInterpolant;
wellInterpolants.roughness = roughnessInterpolant;
% Also build bottom-referenced interpolants F(x_b), x_b in [0,maxDepth], 0 at bottom
xbSampling = maxDepth - depthSampling;            % descending
xbAsc = xbSampling(end:-1:1);                     % ascending from 0->maxDepth
idAsc = innerDiameterSampled(end:-1:1);
thAsc = wallThicknessSampled(end:-1:1);
pfAsc = perforationStatusSampled(end:-1:1);
rgAsc = roughnessSampled(end:-1:1);
wellInterpolants.innerDiameterX = griddedInterpolant(xbAsc, idAsc, 'nearest', 'nearest');
wellInterpolants.wallThicknessX = griddedInterpolant(xbAsc, thAsc, 'nearest', 'nearest');
wellInterpolants.perforationStatusX = griddedInterpolant(xbAsc, pfAsc, 'nearest', 'nearest');
wellInterpolants.roughnessX = griddedInterpolant(xbAsc, rgAsc, 'nearest', 'nearest');
% By default, expose bottom-referenced perforation status under the common name
wellInterpolants.perforationStatus = wellInterpolants.perforationStatusX;
wellInterpolants.offset = offsetInterpolant;
wellInterpolants.angle = angleInterpolant;
wellInterpolants.vertDepth = griddedInterpolant(depthSampling(:), vertDepthSampled(:), 'linear', 'linear');
wellInterpolants.gravityProjection = griddedInterpolant(depthSampling(:), gravProjSampled(:), 'linear', 'nearest');
wellInterpolants.gravityAngle = griddedInterpolant(depthSampling(:), gravAngleSampled(:), 'linear', 'nearest');
% Bottom-referenced deviation/angle
offAsc = offsetSampled(end:-1:1);
angAsc = angleSampled(end:-1:1);
vdAsc = vertDepthSampled(end:-1:1);
gpAsc = gravProjSampled(end:-1:1);
gaAsc = gravAngleSampled(end:-1:1);
wellInterpolants.offsetX = griddedInterpolant(xbAsc, offAsc, 'linear', 'linear');
wellInterpolants.angleX  = griddedInterpolant(xbAsc, angAsc, 'linear', 'linear');
wellInterpolants.vertDepthX = griddedInterpolant(xbAsc, vdAsc, 'linear', 'linear');
wellInterpolants.gravityProjectionX = griddedInterpolant(xbAsc, gpAsc, 'linear', 'nearest');
wellInterpolants.gravityAngleX = griddedInterpolant(xbAsc, gaAsc, 'linear', 'nearest');

% Save to MAT file
save('well_interpolants.mat', 'idInterpolant', 'thicknessInterpolant', 'perfInterpolant', 'roughnessInterpolant', 'offsetInterpolant', 'angleInterpolant');

% Display validation info
%displayValidationInfo(wellInterpolants, minDepth, maxDepth);
end

function displayValidationInfo(interpolants, minDepth, maxDepth)
fprintf('Well Interpolants Created:\n');
fprintf('  Depth range: %.2f m to %.2f m\n', minDepth, maxDepth);

% Verify values at key depths
testDepths = [0, 100, 149, 310, 1000, 1456, 1495, 1634, 1926, 2100];
fprintf('\nValidation at Key Depths:\n');

% Check which interpolants we have
hasDeviationData = ~isempty(interpolants.offset);

if hasDeviationData
    fprintf('Depth (m) | Inner Diam (m) | Wall Thickness (m) | Perforation | Roughness (m) | Offset (m) | Angle (rad)\n');
    fprintf('----------------------------------------------------------------------------------------------------\n');
    for d = testDepths
        id_val = interpolants.innerDiameter(d);
        wt_val = interpolants.wallThickness(d);
        pf_val = interpolants.perforationStatus(d);
        rg_val = interpolants.roughness(d);
        of_val = interpolants.offset(d);
        ang_val = interpolants.angle(d);
        fprintf('%9.2f | %13.4f | %16.4f | %10d | %11.6f | %9.2f | %10.6f\n', d, id_val, wt_val, pf_val, rg_val, of_val, ang_val);
    end
else
    fprintf('Depth (m) | Inner Diam (m) | Wall Thickness (m) | Perforation | Roughness (m)\n');
    fprintf('------------------------------------------------------------------------\n');
    for d = testDepths
        id_val = interpolants.innerDiameter(d);
        wt_val = interpolants.wallThickness(d);
        pf_val = interpolants.perforationStatus(d);
        rg_val = interpolants.roughness(d);
        fprintf('%9.2f | %13.4f | %16.4f | %10d | %11.6f\n', d, id_val, wt_val, pf_val, rg_val);
    end
end
end

function breakpoints = findStepBreakpoints(x, y)
% Find the indices where the function values change
diffY = diff(y);
breakIdx = find(abs(diffY) > eps);
breakpoints = breakIdx;
end

function stepPoints = createStepPoints(x, y, breakpoints)
% Create precise step function representation points
depths = [];
values = [];

% Start with first point
depths = [depths; x(1)];
values = [values; y(1)];

% Add step transitions
for i = 1:length(breakpoints)
    idx = breakpoints(i);

    % Add point just before transition (using same value as before)
    depths = [depths; x(idx+1)-1e-10];
    values = [values; y(idx)];

    % Add point at transition (new value)
    depths = [depths; x(idx+1)];
    values = [values; y(idx+1)];
end

% Add final point
depths = [depths; x(end)];
values = [values; y(end)];

stepPoints = struct('depths', depths, 'values', values);
end

function angleRad = normalizeDeviationAngles(angleValues)
angleRad = angleValues(:);
valid = isfinite(angleRad);
if ~any(valid)
    angleRad(:) = 0;
    return;
end

maxAngle = max(abs(angleRad(valid)));
if maxAngle > pi/2
    fprintf('Converting angle data from degrees to radians\n');
    angleRad(valid) = deg2rad(angleRad(valid));
else
    fprintf('Angle data appears to already be in radians\n');
end

angleRad(~valid) = 0;
end

function [vertDepthSampled, gravProjSampled, gravAngleSampled] = ...
        sampleGravityGeometry(depthSampling, hasDeviationData, deviationData, angleSampled)
depthSampling = depthSampling(:);
angleSampled = angleSampled(:);

vertDepthSampled = depthSampling;
gravProjSampled = ones(size(depthSampling));
gravAngleSampled = zeros(size(depthSampling));

if ~hasDeviationData || isempty(deviationData)
    return;
end

if ismember('VertDepth', deviationData.Properties.VariableNames)
    measDepth = deviationData.MeasDepth(:);
    vertDepth = deviationData.VertDepth(:);
    valid = isfinite(measDepth) & isfinite(vertDepth);
    if nnz(valid) >= 2
        measDepth = measDepth(valid);
        vertDepth = vertDepth(valid);
        [measDepth, uniqIdx] = unique(measDepth, 'stable');
        vertDepth = vertDepth(uniqIdx);
        vertDepthInterpolant = griddedInterpolant(measDepth, vertDepth, 'linear', 'linear');
        vertDepthSampled = vertDepthInterpolant(depthSampling);
    else
        vertDepthSampled = integrateVerticalDepth(depthSampling, angleSampled);
    end
else
    vertDepthSampled = integrateVerticalDepth(depthSampling, angleSampled);
end

gravProjSampled = gradient(vertDepthSampled, depthSampling);
fallbackProj = cos(angleSampled);
fallbackProj(~isfinite(fallbackProj)) = 1;

badProj = ~isfinite(gravProjSampled);
gravProjSampled(badProj) = fallbackProj(badProj);
gravProjSampled = max(-1, min(1, gravProjSampled));
gravAngleSampled = acos(gravProjSampled);
end

function [innerDiameterSampled, wallThicknessSampled, perforationStatusSampled, roughnessSampled] = ...
        sampleCasingProfiles(wellData, depthSampling)
depthSampling = depthSampling(:);
nDepth = numel(depthSampling);

innerDiameterSampled = nan(nDepth, 1);
wallThicknessSampled = nan(nDepth, 1);
roughnessSampled = nan(nDepth, 1);
perforationStatusSampled = false(nDepth, 1);

bestInner = inf(nDepth, 1);
bestOuter = -inf(nDepth, 1);

for i = 1:height(wellData)
    inSection = depthSampling >= wellData.TopMeasDepth(i) & depthSampling <= wellData.BottomMeasDepth(i);
    if ~any(inSection)
        continue;
    end

    innerDiam = wellData.InnerDiameter(i);
    outerDiam = wellData.OuterDiameter(i);

    updateInner = inSection & innerDiam < bestInner;
    innerDiameterSampled(updateInner) = innerDiam;
    roughnessSampled(updateInner) = wellData.Roughness(i);
    bestInner(updateInner) = innerDiam;

    updateOuter = inSection & outerDiam > bestOuter;
    bestOuter(updateOuter) = outerDiam;

    if strcmp(wellData.CasingType{i}, 'Slotted-Liner')
        perforationStatusSampled(inSection) = true;
    end
end

activeDepth = isfinite(bestInner) & isfinite(bestOuter);
wallThicknessSampled(activeDepth) = 0.5 * (bestOuter(activeDepth) - bestInner(activeDepth));
perforationStatusSampled = double(perforationStatusSampled);
end

function vertDepth = integrateVerticalDepth(depthSampling, angleRad)
depthSampling = depthSampling(:);
angleRad = angleRad(:);
vertDepth = zeros(size(depthSampling));

if isempty(depthSampling)
    return;
end

vertDepth(1) = depthSampling(1);
for i = 2:numel(depthSampling)
    ds = depthSampling(i) - depthSampling(i - 1);
    cth = cos(0.5 * (angleRad(i - 1) + angleRad(i)));
    if ~isfinite(cth)
        cth = 1;
    end
    vertDepth(i) = vertDepth(i - 1) + ds * cth;
end
end

% This function has been removed - well plotting is now handled by createWellPlots.m
