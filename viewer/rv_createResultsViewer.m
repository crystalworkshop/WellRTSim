function viewer = rv_createResultsViewer(viewerData)
% rv_createResultsViewer Launch an interactive viewer for simulation results.

    arguments
        viewerData (1, 1) struct
    end

    nTimes = numel(viewerData.times);
    viewer = struct();
    viewer.data = viewerData;
    viewer.cache = cell(1, nTimes);

    fig = uifigure('Name', 'Simulation Results Viewer', ...
        'Color', [0.12 0.12 0.14], ...
        'AutoResizeChildren', 'on', ...
        'Resize', 'on');

    screenSize = get(0, 'ScreenSize');
    figWidth = min(1200, max(920, round(0.78 * screenSize(3))));
    figHeight = min(760, max(540, round(0.68 * screenSize(4))));
    fig.Position(3:4) = [figWidth, figHeight];
    fig.Position(1:2) = [(screenSize(3) - figWidth) / 2, (screenSize(4) - figHeight) / 2];

    mainLayout = uigridlayout(fig, [2, 1]);
    mainLayout.RowHeight = {'1x', 80};
    mainLayout.Padding = [6 6 6 6];
    mainLayout.BackgroundColor = fig.Color;

    axesContainer = uipanel(mainLayout, 'BorderType', 'none');
    axesContainer.Layout.Row = 1;
    axesContainer.BackgroundColor = fig.Color;
    axesStruct = rv_setupResultsViewerAxes(axesContainer);

    exportButton = [];
    if isfield(axesStruct, 'exportScaleButton')
        exportButton = axesStruct.exportScaleButton;
        exportButton.ButtonPushedFcn = @(~, ~) exportScaleFractions();
    end

    ctrlLayout = uigridlayout(mainLayout, [2, 3]);
    ctrlLayout.Layout.Row = 2;
    ctrlLayout.ColumnWidth = {'fit', '1x', 'fit'};
    ctrlLayout.RowHeight = {'fit', 'fit'};
    ctrlLayout.Padding = [10 5 10 5];
    ctrlLayout.BackgroundColor = fig.Color;

    fileLabel = uilabel(ctrlLayout, ...
        'Text', sprintf('File: %s', viewerData.file), ...
        'FontColor', [0.85 0.85 0.85], ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'WordWrap', 'on');
    fileLabel.Layout.Row = 1;
    fileLabel.Layout.Column = [1 3];

    sliderLabel = uilabel(ctrlLayout, ...
        'Text', sprintf('Profile time (%s)', viewerData.timeUnitLabel), ...
        'FontColor', [0.85 0.85 0.85], ...
        'FontSize', 10);
    sliderLabel.Layout.Row = 2;
    sliderLabel.Layout.Column = 1;

    sliderMax = max(1, nTimes);
    sliderLimits = [1, max(2, sliderMax)];
    slider = uislider(ctrlLayout, ...
        'Limits', sliderLimits, ...
        'Value', sliderMax, ...
        'ValueChangingFcn', @(~, evt) previewIndex(evt.Value), ...
        'ValueChangedFcn', @(src, ~) commitIndex(src.Value));
    slider.Layout.Row = 2;
    slider.Layout.Column = 2;
    slider.MajorTicks = computeTickPositions(nTimes);
    slider.MajorTickLabels = formatTickLabels(slider.MajorTicks, viewerData);
    if nTimes <= 1
        slider.Enable = 'off';
    end

    timeLabel = uilabel(ctrlLayout, ...
        'Text', '', ...
        'FontColor', [0.9 0.9 0.9], ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'right');
    timeLabel.Layout.Row = 2;
    timeLabel.Layout.Column = 3;

    viewer.fig = fig;
    viewer.axes = axesStruct;
    viewer.slider = slider;
    viewer.timeLabel = timeLabel;
    viewer.depthSlider = axesStruct.depthSlider;
    viewer.previewIndex = max(1, nTimes);
    viewer.currentTimeIndex = max(1, nTimes);
    viewer.currentDepthIndex = max(1, min(numel(viewerData.x), viewerData.defaultDepthIndex));
    viewer.hasScaleData = any(viewerData.scaleThickness(:) > 0);

    if ~isempty(exportButton)
        if isempty(viewerData.scaleFractions) || ~viewer.hasScaleData
            exportButton.Enable = 'off';
            exportButton.Tooltip = 'No scale fraction data available';
        else
            exportButton.Enable = 'on';
        end
    end

    configureDepthSlider();

    if nTimes > 0
        drawProfile(nTimes);
    else
        timeLabel.Text = 'No time steps available';
    end

    function drawProfile(idx)
        idx = max(1, min(numel(viewerData.times), round(idx)));
        viewer.previewIndex = idx;
        if isempty(viewer.cache{idx})
            viewer.cache{idx} = rv_readProfileSnapshot(viewerData, idx);
        end
        profile = viewer.cache{idx};
        viewer.currentTimeIndex = idx;
        viewer.currentProfile = profile;
        rv_updateResultsViewerPlots(axesStruct, viewerData, profile, idx, viewer.currentDepthIndex);
        displayTime = profile.timeDays * viewerData.timeScale;
        timeLabel.Text = sprintf('t = %.3g %s', displayTime, viewerData.timeUnitLabel);
    end

    function previewIndex(val)
        if isempty(viewerData.times)
            return;
        end
        drawProfile(val);
    end

    function commitIndex(val)
        drawProfile(val);
        slider.Value = viewer.previewIndex;
    end

    function configureDepthSlider()
        depthSlider = viewer.depthSlider;
        if ~viewer.hasScaleData
            depthSlider.Enable = 'off';
            return;
        end

        nDepth = numel(viewerData.x);
        minDepth = min(viewerData.x);
        maxDepth = max(viewerData.x);
        depthSlider.Enable = 'on';
        depthSlider.Limits = [minDepth maxDepth];
        tickDepths = linspace(minDepth, maxDepth, min(6, nDepth));
        tickDepths = unique(tickDepths, 'stable');
        depthSlider.MajorTicks = tickDepths;
        depthSlider.MajorTickLabels = compose('%.0f', tickDepths);
        viewer.currentDepthIndex = max(1, min(nDepth, viewer.currentDepthIndex));
        depthSlider.Value = viewerData.x(viewer.currentDepthIndex);
        depthSlider.ValueChangingFcn = @(~, evt) depthSliderChanging(evt.Value);
        depthSlider.ValueChangedFcn = @(src, ~) depthSliderChanged(src.Value);
        depthSlider.Tooltip = 'Select depth [m]';
    end

    function depthSliderChanging(val)
        refreshScaleDepth(val, false);
    end

    function depthSliderChanged(val)
        refreshScaleDepth(val, true);
    end

    function refreshScaleDepth(val, commit)
        if ~viewer.hasScaleData
            return;
        end
        targetDepth = max(min(viewerData.x), min(max(viewerData.x), val));
        if commit
            viewer.depthSlider.Value = targetDepth;
        end
        [~, idx] = min(abs(viewerData.x - targetDepth));
        viewer.currentDepthIndex = idx;

        if viewer.currentTimeIndex < 1 || viewer.currentTimeIndex > numel(viewerData.times)
            return;
        end
        profile = viewer.cache{viewer.currentTimeIndex};
        if isempty(profile)
            profile = rv_readProfileSnapshot(viewerData, viewer.currentTimeIndex);
            viewer.cache{viewer.currentTimeIndex} = profile;
        end
        rv_updateResultsViewerPlots(axesStruct, viewerData, profile, ...
            viewer.currentTimeIndex, viewer.currentDepthIndex);
    end

    function exportScaleFractions()
        if isempty(viewerData.scaleFractions)
            uialert(fig, 'No scale fraction data available.', 'Export Scale Fractions');
            return;
        end

        nMinerals = size(viewerData.scaleFractions, 1);
        nSeriesTimes = numel(viewerData.times);
        if nMinerals == 0 || nSeriesTimes == 0
            uialert(fig, 'No scale fraction data available.', 'Export Scale Fractions');
            return;
        end

        depthVal = viewer.depthSlider.Value;
        if isempty(depthVal) || ~isfinite(depthVal)
            depthIdx = max(1, min(numel(viewerData.x), viewer.currentDepthIndex));
        else
            [~, depthIdx] = min(abs(viewerData.x - depthVal));
        end
        depthVal = viewerData.x(depthIdx);

        fracSeries = reshape(viewerData.scaleFractions(:, depthIdx, :), nMinerals, nSeriesTimes);
        thicknessSeries = viewerData.scaleThickness(depthIdx, :).';
        timeVals = viewerData.times(:) * viewerData.timeScale;
        timeVar = matlab.lang.makeValidName(sprintf('time_%s', char(viewerData.timeUnitLabel)));
        T = table(timeVals, thicknessSeries, 'VariableNames', {timeVar, 'scaleThickness_m'});

        names = getScaleFractionNames(nMinerals, viewerData);
        for k = 1:nMinerals
            T.(matlab.lang.makeValidName(['scaleFrac_' names{k}])) = fracSeries(k, :).';
        end

        outDir = fileparts(viewerData.file);
        outFile = fullfile(outDir, sprintf('scale_%s.csv', formatDepthToken(depthVal)));
        try
            writetable(T, outFile);
            uialert(fig, sprintf('Saved scale fractions at depth %.1f m to:\n%s', ...
                depthVal, outFile), 'Scale Fractions Saved');
        catch ME
            uialert(fig, sprintf('Failed to write %s\n%s', outFile, ME.message), ...
                'Export Failed');
        end
    end
end

function token = formatDepthToken(depthVal)
    token = sprintf('%.1f', depthVal);
    token = regexprep(token, '\.?0+$', '');
end

function names = getScaleFractionNames(nMinerals, viewerData)
    names = cell(1, nMinerals);
    raw = {};
    if isfield(viewerData, 'scaleMineralNames') && ~isempty(viewerData.scaleMineralNames)
        raw = cellstr(viewerData.scaleMineralNames(:));
    elseif isfield(viewerData, 'scaleFractionDatasets') && ~isempty(viewerData.scaleFractionDatasets)
        raw = viewerData.scaleFractionDatasets(:);
    end
    for k = 1:nMinerals
        if k <= numel(raw) && ~isempty(raw{k})
            names{k} = char(raw{k});
        else
            names{k} = sprintf('Mineral%d', k);
        end
    end
end

function ticks = computeTickPositions(n)
    if n <= 5
        ticks = 1:n;
    else
        step = max(1, round(n / 5));
        ticks = unique([1, step:step:n, n]);
    end
end

function labels = formatTickLabels(indices, viewerData)
    if isempty(indices)
        labels = {};
        return;
    end
    indices = round(indices);
    indices(indices < 1) = 1;
    indices(indices > numel(viewerData.times)) = numel(viewerData.times);
    labels = compose('%.3g', viewerData.times(indices) * viewerData.timeScale);
end
