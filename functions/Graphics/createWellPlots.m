function createWellPlots(state, wellData, deviationData, minDepth, maxDepth, wellName)
% Create well schematic and trajectory plots in the designated tab
%
% Inputs:
%   state - State structure with graphics components
%   wellData - Well casing data table
%   deviationData - Well deviation data table  
%   minDepth - Minimum depth for plotting
%   maxDepth - Maximum depth for plotting
%   wellName - Name of the well

% Ensure tab is ready for content
drawnow;

% Rebuild the tab from a clean slate when graphics are reinitialized.
delete(state.schematicTab.Children);

% Create grid layout with 2 equal columns that fills the entire tab
gl = uigridlayout(state.schematicTab, [1, 2]); % 1 rows 2 columns for balanced layout
gl.ColumnWidth = {'1x', '1x'}; % Equal column widths
gl.RowHeight = {'1x'}; % Use full height  
gl.Padding = [10, 10, 10, 10];
gl.ColumnSpacing = 15; % Space between columns

% Grid layout will automatically fill the tab - no explicit positioning needed

% Create two axes using the grid layout with explicit positioning
axLeft = uiaxes(gl);
axLeft.Layout.Row = 1;
axLeft.Layout.Column = 1;

axRight = uiaxes(gl);
axRight.Layout.Row = 1;
axRight.Layout.Column = 2;

% Force MATLAB to update the display and layout before plotting
drawnow;

% Force the grid layout to re-calculate its layout
gl.ColumnWidth = {'1x', '1x'}; % Reset column widths to trigger recalculation
drawnow; % Additional drawnow to ensure layout is applied

% === LEFT PANEL: WELL SCHEMATIC ===
% Define plot area for vertical schematic (centered)
xMin = -0.3;  % Center the well plot
xMax = 0.3;
yMin = minDepth;
yMax = maxDepth;

% Set up axes properties and limits first
hold(axLeft, 'on');
hold(axRight, 'on');
axis(axLeft, [xMin xMax yMin yMax]);
set(axLeft, 'YDir', 'reverse');  % Depth increases downward

% Set up right axis limits early
maxOffset = 100;  % Default width
ylim(axRight, [yMin yMax]);
xlim(axRight, [-maxOffset, maxOffset]);
set(axRight, 'YDir', 'reverse');
title(axLeft, [wellName, ' Geothermal Well Schematic']);
xlabel(axLeft, 'Diameter (m)');
ylabel(axLeft, 'Depth (m)');
grid(axLeft, 'on');

% Plot surface (across the full width)
plot(axLeft, [xMin, xMax], [0, 0], 'k-', 'LineWidth', 2);

% Sort the wellData by diameter (largest first) to draw from outside in
wellData = sortrows(wellData, 'OuterDiameter', 'descend');

% Pattern for slotted liners
slotPattern = 20;  % Spacing for slots in points

% First draw the outer casings as shapes
for i = 1:height(wellData)
    casing = wellData(i, :);
    topDepth = casing.TopMeasDepth;
    bottomDepth = casing.BottomMeasDepth;
    outerDiam = casing.OuterDiameter;
    innerDiam = casing.InnerDiameter;
    casingType = casing.CasingType{1};

    % Plot outer casing outline
    if strcmp(casingType, 'Slotted-Liner')
        % For slotted liner, draw with pattern
        y = linspace(topDepth, bottomDepth, ceil((bottomDepth-topDepth)/slotPattern));
        for j = 1:length(y)-1
            % Left side slots
            plot(axLeft, [-outerDiam/2, -innerDiam/2], [y(j), y(j)], 'k-', 'LineWidth', 1);
            % Right side slots
            plot(axLeft, [outerDiam/2, innerDiam/2], [y(j), y(j)], 'k-', 'LineWidth', 1);
        end
    end

    % Draw casing outline
    x = [-outerDiam/2, outerDiam/2, outerDiam/2, -outerDiam/2, -outerDiam/2];
    y = [topDepth, topDepth, bottomDepth, bottomDepth, topDepth];

    % Fill with color based on casing type
    fill(axLeft, x, y, getCasingColor(casingType), 'EdgeColor', 'k', 'LineWidth', 1);

    % Draw inner hole
    x = [-innerDiam/2, innerDiam/2, innerDiam/2, -innerDiam/2, -innerDiam/2];
    fill(axLeft, x, y, 'w', 'EdgeColor', 'k', 'LineWidth', 1);
end

% Get all unique transition depths from the data
depthMarkers = unique([wellData.TopMeasDepth; wellData.BottomMeasDepth]);
depthMarkers = sort(depthMarkers);

% Add depth markers in the center of the wellbore
for i = 1:length(depthMarkers)
    d = depthMarkers(i);
    if d >= minDepth && d <= maxDepth
        % Place text in the center of the wellbore (x=0) with white background
        text(axLeft, 0, d, sprintf('%.0f m', d), ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'Color', 'black', ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'black', ...
            'LineWidth', 0.5, ...
            'Margin', 3);
    end
end

% Add legend for vertical schematic
h1 = plot(axLeft, NaN, NaN, 'Color', getCasingColor('Casing'), 'LineWidth', 10);
h2 = plot(axLeft, NaN, NaN, 'Color', getCasingColor('Liner'), 'LineWidth', 10);
h3 = plot(axLeft, NaN, NaN, 'Color', getCasingColor('Slotted-Liner'), 'LineWidth', 10);
legend(axLeft, [h1, h2, h3], {'Casing', 'Liner', 'Slotted-Liner'}, 'Location', 'best');

% Add annotations for casing sizes dynamically
% Calculate label positions for each casing type (approx 1/3 down from top)
for i = 1:height(wellData)
    casing = wellData(i, :);
    labelDepth = casing.TopMeasDepth + (casing.BottomMeasDepth - casing.TopMeasDepth) / 3;
    labelText = sprintf('%s %s', casing.CasingSize{1}, casing.CasingType{1});
    text(axLeft, xMax-0.05, labelDepth, labelText, 'HorizontalAlignment', 'right');
end

% === RIGHT PANEL: WELL TRAJECTORY ===
% Create a smooth trajectory from the deviation data
[measDepths, offsets] = computeTrajectoryOffsets(deviationData);

% Create interpolants for smooth plotting
offsetInterp = griddedInterpolant(measDepths, offsets, 'linear', 'linear');

% Sample points for smooth curve
depthSamples = linspace(min(measDepths), max(measDepths), 1000);
offsetSamples = offsetInterp(depthSamples);

% Plot well trajectory
hold(axRight, 'on');
hPath = plot(axRight, offsetSamples, depthSamples, 'b-', 'LineWidth', 2);

% Add surface
maxOffset = max(abs(offsetSamples)) * 1.1;
if maxOffset < 10
    maxOffset = 100;  % Minimum display width for clarity
end

% Format trajectory plot
grid(axRight, 'on');
title(axRight, [wellName, ' Well Trajectory']);
xlabel(axRight, 'Horizontal Offset (m)');
ylabel(axRight, 'Depth (m)');
set(axRight, 'YDir', 'reverse');
ylim(axRight, [minDepth, maxDepth]);
xlim(axRight, [-maxOffset, maxOffset]);
axis(axRight, 'equal');

% Mark slotted liner sections on the trajectory
hasSlottedLiner = false;
hSlotted = [];
for i = 1:height(wellData)
    casing = wellData(i, :);
    if strcmp(casing.CasingType{1}, 'Slotted-Liner')
        hasSlottedLiner = true;
        % Get a few points along the slotted liner segment
        slottedDepths = linspace(casing.TopMeasDepth, casing.BottomMeasDepth, 10);
        slottedOffsets = offsetInterp(slottedDepths);

        % Draw markers for slotted liner
        h = plot(axRight, slottedOffsets, slottedDepths, 'c-', 'LineWidth', 4);
        if isempty(hSlotted)
            hSlotted = h;
        end

        % Label the slotted section
        midDepth = (casing.TopMeasDepth + casing.BottomMeasDepth) / 2;
        midOffset = offsetInterp(midDepth);
        text(axRight, midOffset + maxOffset*0.05, midDepth, 'Slotted Liner', 'Color', 'blue');
    end
end

% Add feed zone intervals along the trajectory if available
hFeed = [];
if isfield(state, 'feedzones')
    fz = state.feedzones;
    depthMin = [];
    depthMax = [];
    if isfield(fz, 'depth_min') && isfield(fz, 'depth_max') && ~isempty(fz.depth_min)
        depthMin = fz.depth_min(:);
        depthMax = fz.depth_max(:);
    elseif isfield(fz, 'depth') && ~isempty(fz.depth)
        midDepths = fz.depth(:);
        segHalf = 0.5;
        if isfield(state, 'dx') && isfinite(state.dx) && state.dx > 0
            segHalf = max(segHalf, 0.5 * state.dx);
        end
        depthMin = midDepths - segHalf;
        depthMax = midDepths + segHalf;
    end

    if ~isempty(depthMin)
        segHalf = 0.5;
        if isfield(state, 'dx') && isfinite(state.dx) && state.dx > 0
            segHalf = max(segHalf, 0.5 * state.dx);
        end
        for k = 1:numel(depthMin)
            d1 = depthMin(k);
            d2 = depthMax(k);
            if ~isfinite(d1) || ~isfinite(d2)
                continue;
            end
            dStart = max(minDepth, min(d1, d2));
            dEnd = min(maxDepth, max(d1, d2));
            if dEnd <= dStart
                dStart = max(minDepth, dStart - segHalf);
                dEnd = min(maxDepth, dEnd + segHalf);
            end
            if dEnd <= dStart
                continue;
            end
            nSeg = max(2, ceil(abs(dEnd - dStart) / 5));
            depthSeg = linspace(dStart, dEnd, nSeg);
            offsetSeg = offsetInterp(depthSeg);
            h = plot(axRight, offsetSeg, depthSeg, '-', 'LineWidth', 4, 'Color', [1.0 0.9 0.1]);
            if isempty(hFeed)
                hFeed = h;
            end
        end
    end
end

% Add a legend for the trajectory with optimal placement
handles = hPath;
labels = {'Wellbore Path'};
if ~isempty(hSlotted)
    handles(end+1) = hSlotted; %#ok<AGROW>
    labels{end+1} = 'Production Zone'; %#ok<AGROW>
end
if ~isempty(hFeed)
    handles(end+1) = hFeed; %#ok<AGROW>
    labels{end+1} = 'Feedzone Intervals'; %#ok<AGROW>
end
leg = legend(axRight, handles, labels, 'Location', 'eastoutside');
leg.FontSize = 9;
leg.Box = 'on';

% Final update to ensure proper display
drawnow;

% Force axes to properly fit their content
axis(axLeft, 'tight');
ylim(axLeft, [yMin yMax]);
xlim(axLeft, [xMin xMax]);

axis(axRight, 'tight');
ylim(axRight, [yMin yMax]); 

% Final drawnow to ensure everything is rendered
drawnow;

end

function color = getCasingColor(casingType)
switch strtrim(casingType)
    case 'Casing'
        color = [0.7, 0.7, 0.7];
    case 'Liner'
        color = [0.5, 0.5, 0.5];
    case 'Slotted-Liner'
        color = [0.3, 0.6, 0.9];
    otherwise
        error('createWellPlots:UnknownCasingType', ...
            'Unknown casing type "%s" in well plot data.', casingType);
end
end

function [measDepths, offsets] = computeTrajectoryOffsets(deviationData)
% Compute offsets from angle + vertical depth, accumulating from the top.
measDepths = deviationData.MeasDepth(:);
offsetsRaw = deviationData.Offset(:);
vertDepth = measDepths;
if ismember('VertDepth', deviationData.Properties.VariableNames)
    vertDepth = deviationData.VertDepth(:);
end
angle = deviationData.Angle(:);

n = min([numel(measDepths), numel(offsetsRaw), numel(vertDepth), numel(angle)]);
measDepths = measDepths(1:n);
offsetsRaw = offsetsRaw(1:n);
vertDepth = vertDepth(1:n);
angle = angle(1:n);

if isempty(measDepths)
    offsets = offsetsRaw;
    return;
end

if measDepths(1) > 0
    measDepths = [0; measDepths];
    vertDepth = [0; vertDepth];
    angle = [angle(1); angle];
    offsetsRaw = [0; offsetsRaw];
end

validAngles = angle(isfinite(angle));
if isempty(validAngles)
    offsets = zeros(size(measDepths));
    return;
end

maxAngle = max(abs(validAngles));
if maxAngle > pi/2
    angleRad = deg2rad(angle);
else
    angleRad = angle;
end

offsets = offsetsRaw;
if ~isfinite(offsets(1))
    offsets(1) = 0;
end

for i = 2:numel(measDepths)
    if isfinite(offsetsRaw(i)) && offsetsRaw(i) ~= 0
        offsets(i) = offsetsRaw(i);
        continue;
    end
    dvert = vertDepth(i) - vertDepth(i-1);
    ang = angleRad(i);
    if ~isfinite(dvert) || ~isfinite(ang)
        dvert = 0;
        ang = 0;
    end
    offsets(i) = offsets(i-1) + dvert * tan(ang);
end
end
