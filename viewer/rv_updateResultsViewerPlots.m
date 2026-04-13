function rv_updateResultsViewerPlots(axesStruct, viewerData, profile, timeIndex, depthIndex)
% rv_updateResultsViewerPlots Render one profile snapshot on prepared axes.

    arguments
        axesStruct (1, 1) struct
        viewerData (1, 1) struct
        profile (1, 1) struct
        timeIndex (1, 1) double {mustBePositive, mustBeInteger} = 1
        depthIndex (1, 1) double {mustBePositive, mustBeInteger} = 1
    end

    x = viewerData.x(:);
    if numel(x) ~= numel(profile.P)
        x = x(1:min(numel(x), numel(profile.P)));
    end

    displayTime = profile.timeDays * viewerData.timeScale;
    timeStr = sprintf('%.3g %s', displayTime, viewerData.timeUnitLabel);

    ax = axesStruct.ax1;
    cla(ax);
    hold(ax, 'on');
    [xP, pData] = alignVectors(x, profile.P);
    plot(ax, xP, pData, '-', 'LineWidth', 1.2, 'Color', [0.0 1.0 1.0]);
    if hasMeasured(viewerData, 'measuredPressure')
        plotMeasuredProfile(ax, viewerData.measuredPressure, [0.95 0.8 0.1], 'MeasuredPressure');
    end
    ax.Title.String = sprintf('Pressure  t = %s', timeStr);
    ax.XLabel.String = 'Position [m]';
    ax.YLabel.String = sprintf('P [%s]', viewerData.pressureUnitLabel);
    hold(ax, 'off');

    ax = axesStruct.ax2;
    cla(ax);
    hold(ax, 'on');
    [xT, tData] = alignVectors(x, profile.T_C);
    plot(ax, xT, tData, '-', 'LineWidth', 1.2, 'Color', [1.0 0.2 0.2]);
    if hasMeasured(viewerData, 'measuredTemperature')
        plotMeasuredProfile(ax, viewerData.measuredTemperature, [0.4 0.8 1.0], 'MeasuredTemperature');
    end
    ax.Title.String = sprintf('Temperature  t = %s', timeStr);
    ax.XLabel.String = 'Position [m]';
    ax.YLabel.String = 'T [degC]';
    hold(ax, 'off');

    ax = axesStruct.ax3;
    cla(ax);
    hold(ax, 'on');
    [xUmix, umix] = alignVectors(x, profile.u_mix);
    [xUgas, ugas] = alignVectors(x, profile.u_gas);
    [xUliq, uliq] = alignVectors(x, profile.u_liq);
    legendHandles = [];
    legendLabels = {};
    if ~isempty(umix)
        h = plot(ax, xUmix, umix, '-', 'LineWidth', 1.2, 'Color', [0.2 0.8 0.2]);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'Mixture'; %#ok<AGROW>
    end
    if ~isempty(ugas)
        h = plot(ax, xUgas, ugas, '--', 'LineWidth', 1.0, 'Color', [1.0 0.5 0.5]);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'Gas'; %#ok<AGROW>
    end
    if ~isempty(uliq)
        h = plot(ax, xUliq, uliq, '--', 'LineWidth', 1.0, 'Color', [0.4 0.9 0.9]);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'Liquid'; %#ok<AGROW>
    end
    if hasMeasuredSeries(viewerData, 'measuredVelocity')
        [mh, mlabels] = plotMeasuredSeries(ax, viewerData.measuredVelocity, ...
            {[0.95 0.7 0.1], [0.65 0.6 1.0]}, 'MeasuredVel');
        legendHandles = [legendHandles, mh]; %#ok<AGROW>
        legendLabels = [legendLabels, mlabels]; %#ok<AGROW>
    end
    if ~isempty(legendHandles)
        leg = legend(ax, legendHandles, legendLabels, ...
            'Location', 'southoutside', 'Orientation', 'horizontal');
        applyLegendStyle(leg);
    else
        legend(ax, 'off');
    end
    ax.Title.String = sprintf('Velocities  t = %s', timeStr);
    ax.XLabel.String = 'Position [m]';
    ax.YLabel.String = 'u [m/s]';
    hold(ax, 'off');

    renderFlowRates(axesStruct.ax6, viewerData, profile);

    ax = axesStruct.ax4;
    cla(ax);
    hold(ax, 'on');
    plot(ax, x, zeros(size(x)), '--', 'LineWidth', 1.0, 'Color', [0.6 0.6 0.6], ...
        'HandleVisibility', 'off');
    if ~isempty(profile.SI)
        colors = lines(max(numel(profile.SI), 1));
        legendEntries = {};
        legendHandles = [];
        for k = 1:numel(profile.SI)
            data = profile.SI{k};
            if isempty(data)
                continue;
            end
            [xSi, siData] = alignVectors(x, data);
            h = plot(ax, xSi, siData, '-', 'LineWidth', 1.1, 'Color', colors(k, :));
            legendHandles(end+1) = h; %#ok<AGROW>
            if k <= numel(viewerData.siNames)
                legendEntries{end+1} = char(viewerData.siNames(k)); %#ok<AGROW>
            else
                legendEntries{end+1} = sprintf('SI %d', k); %#ok<AGROW>
            end
        end
        if ~isempty(legendHandles)
            leg = legend(ax, legendHandles, legendEntries, ...
                'Location', 'southoutside', 'Orientation', 'horizontal', ...
                'NumColumns', max(1, ceil(numel(legendEntries) / 2)));
            applyLegendStyle(leg);
        end
    else
        text(ax, 0.5, 0.5, 'No saturation data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
    end
    ax.Title.String = sprintf('Saturation Index  t = %s', timeStr);
    ax.XLabel.String = 'Position [m]';
    ax.YLabel.String = 'SI [-]';
    ax.YLim = softYLimitLinear(profile.SI, [-2 2]);
    hold(ax, 'off');

    ax = axesStruct.ax5;
    cla(ax);
    ax.YScale = 'log';
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    hold(ax, 'on');
    ppmProfiles = zeros(0, 0);
    if ~isempty(profile.elements)
        [ppmProfiles, xChem] = computeElementPpmProfiles(x, profile.elements);
        colors = lines(max(size(ppmProfiles, 1), 1));
        legendEntries = {};
        legendHandles = [];
        for k = 1:size(ppmProfiles, 1)
            data = ppmProfiles(k, :);
            if isempty(data) || ~any(isfinite(data))
                continue;
            end
            h = semilogy(ax, xChem, data, '-', 'LineWidth', 1.1, 'Color', colors(k, :));
            legendHandles(end+1) = h; %#ok<AGROW>
            if k <= numel(viewerData.elementNames)
                legendEntries{end+1} = char(viewerData.elementNames(k)); %#ok<AGROW>
            else
                legendEntries{end+1} = sprintf('Elem %d', k); %#ok<AGROW>
            end
        end
        if ~isempty(legendHandles)
            leg = legend(ax, legendHandles, legendEntries, ...
                'Location', 'southoutside', 'Orientation', 'horizontal', ...
                'NumColumns', max(1, ceil(numel(legendEntries) / 3)));
            applyLegendStyle(leg);
        end
    else
        text(ax, 0.5, 0.5, 'No concentration data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
    end
    ax.YLim = softYLimitLog(ppmProfiles, [1e-2 1e5]);
    ax.Title.String = sprintf('Mass Concentration  t = %s', timeStr);
    ax.XLabel.String = 'Position [m]';
    ax.YLabel.String = 'c [ppm]';
    hold(ax, 'off');

    ax = axesStruct.axGeom;
    cla(ax);
    hold(ax, 'on');
    [currGeomX, currGeomDp] = alignVectors(x, profile.Dp);
    [baseGeomX, baseGeomDp] = alignVectors(x, viewerData.Dp0);
    if ~isempty(currGeomDp) && ~isempty(baseGeomDp)
        geomCount = min(numel(currGeomDp), numel(baseGeomDp));
        geomPos = currGeomX(1:geomCount);
        currSlice = currGeomDp(1:geomCount);
        baseSlice = baseGeomDp(1:geomCount);
        palette = getMineralPalette(viewerData);
        fracSliceFull = computeScaleFractionSlice(viewerData, timeIndex);
        if isempty(fracSliceFull)
            fracSlice = [];
        elseif size(fracSliceFull, 2) < geomCount
            fracSlice = zeros(size(fracSliceFull, 1), geomCount);
            fracSlice(:, 1:size(fracSliceFull, 2)) = fracSliceFull;
        else
            fracSlice = fracSliceFull(:, 1:geomCount);
        end
        plotGeometryScaleFill(ax, geomPos, baseSlice, currSlice, fracSlice, palette);
    end
    if ~isempty(currGeomDp)
        plot(ax, currGeomDp, currGeomX, '-', 'LineWidth', 1.2, 'Color', [0.2 0.8 0.2]);
    end
    if ~isempty(baseGeomDp)
        plot(ax, baseGeomDp, baseGeomX, '--', 'LineWidth', 1.0, 'Color', [0.7 0.7 0.7]);
    end
    ax.Title.String = sprintf('Geometry  t = %s', timeStr);
    ax.XLabel.String = 'Diameter [m]';
    ax.YLabel.String = 'Position [m]';
    if ~isempty(x)
        ax.YLim = [min(x) max(x)];
    end
    if ~isempty(viewerData.x)
        idx = max(1, min(numel(viewerData.x), round(depthIndex)));
        depthVal = viewerData.x(idx);
        depthLine = getappdata(ax, 'ScaleDepthLine');
        if isempty(depthLine) || ~isvalid(depthLine)
            depthLine = yline(ax, depthVal, '--', ...
                'Color', [1.0 0.6 0.0], ...
                'LineWidth', 1.2, ...
                'Tag', 'ScaleDepthLine', ...
                'HandleVisibility', 'off');
            if isprop(depthLine, 'HitTest')
                depthLine.HitTest = 'off';
            end
            setappdata(ax, 'ScaleDepthLine', depthLine);
        else
            depthLine.Value = depthVal;
            depthLine.Visible = 'on';
        end
    end
    hold(ax, 'off');

    if isvalid(axesStruct.axScale)
        renderScaleCompositionPlot(axesStruct.axScale, viewerData, depthIndex, timeIndex);
    end
    if isvalid(axesStruct.axScaleThickness)
        renderScaleThicknessHistory(axesStruct.axScaleThickness, viewerData, depthIndex, timeIndex);
    end
    if isvalid(axesStruct.axFlow)
        renderPrecipitationRate(axesStruct.axFlow, viewerData, timeIndex, depthIndex);
    end
end

function [ppmProfiles, xChem] = computeElementPpmProfiles(x, elementCells)
    ppmProfiles = zeros(0, 0);
    xChem = [];
    if isempty(elementCells)
        return;
    end

    nElem = numel(elementCells);
    lengths = zeros(1, nElem);
    for k = 1:nElem
        lengths(k) = numel(elementCells{k});
    end
    minLen = min([numel(x), lengths(lengths > 0)]);
    if isempty(minLen) || minLen == 0
        return;
    end

    xChem = x(1:minLen);
    conc = nan(nElem, minLen);
    for k = 1:nElem
        data = elementCells{k};
        if isempty(data)
            continue;
        end
        conc(k, :) = data(1:minLen);
    end

    waterMassfrac = 1 - sum(conc, 1, 'omitnan');
    ppmProfiles = 1e6 * conc ./ waterMassfrac;
    ppmProfiles(:, waterMassfrac <= 0) = NaN;
end

function ylimOut = softYLimitLinear(data, softBounds)
    values = collectFiniteValues(data);
    if isempty(values)
        ylimOut = softBounds;
        return;
    end

    dataMin = min(values);
    dataMax = max(values);
    if dataMin < softBounds(1) || dataMax > softBounds(2)
        ylimOut = softBounds;
        return;
    end

    span = dataMax - dataMin;
    if span < 1e-9
        pad = max(0.1, 0.1 * max(1, abs(dataMax)));
    else
        pad = 0.08 * span;
    end
    ylimOut = [dataMin - pad, dataMax + pad];
    ylimOut(1) = max(ylimOut(1), softBounds(1));
    ylimOut(2) = min(ylimOut(2), softBounds(2));
    if ylimOut(1) >= ylimOut(2)
        ylimOut = softBounds;
    end
end

function ylimOut = softYLimitLog(data, softBounds)
    values = collectFiniteValues(data);
    values = values(values > 0);
    if isempty(values)
        ylimOut = softBounds;
        return;
    end

    dataMin = min(values);
    dataMax = max(values);
    if dataMin < softBounds(1) || dataMax > softBounds(2)
        ylimOut = softBounds;
        return;
    end

    logMin = log10(dataMin);
    logMax = log10(dataMax);
    if logMax - logMin < 1e-6
        logPad = 0.25;
    else
        logPad = 0.08 * (logMax - logMin);
    end

    ylimOut = [10^(logMin - logPad), 10^(logMax + logPad)];
    ylimOut(1) = max(ylimOut(1), softBounds(1));
    ylimOut(2) = min(ylimOut(2), softBounds(2));
    if ylimOut(1) >= ylimOut(2)
        ylimOut = softBounds;
    end
end

function values = collectFiniteValues(data)
    if iscell(data)
        values = zeros(0, 1);
        for i = 1:numel(data)
            values = [values; collectFiniteValues(data{i})]; %#ok<AGROW>
        end
        return;
    end

    if isempty(data)
        values = zeros(0, 1);
        return;
    end

    values = data(:);
    values = values(isfinite(values));
end

function [xOut, qOut] = computePhaseFlow(x, profile, phase)
    xOut = [];
    qOut = [];
    switch string(lower(phase))
        case "gas"
            if ~isempty(profile.Q_v)
                [xOut, qOut] = alignVectors(x, profile.Q_v);
                return;
            end
            if ~isempty(profile.Qm_v)
                [xOut, qOut] = alignVectors(x, profile.Qm_v);
                return;
            end
        otherwise
            if ~isempty(profile.Q_l)
                [xOut, qOut] = alignVectors(x, profile.Q_l);
                return;
            end
            if ~isempty(profile.Qm_l)
                [xOut, qOut] = alignVectors(x, profile.Qm_l);
                return;
            end
    end
end

function tf = hasMeasured(viewerData, fieldName)
    tf = false;
    data = viewerData.(fieldName);
    tf = isstruct(data) && ~isempty(data.x) && ~isempty(data.value);
end

function renderFlowRates(ax, viewerData, profile)
    cla(ax);
    hold(ax, 'on');

    x = viewerData.x(:);
    legendHandles = [];
    legendLabels = {};

    [xGas, qGas] = computePhaseFlow(x, profile, "gas");
    [xLiq, qLiq] = computePhaseFlow(x, profile, "liquid");
    [xMix, qMix] = sumPhaseFlows(profile, x, xGas, qGas, xLiq, qLiq);

    if ~isempty(qMix)
        h = plot(ax, xMix, qMix, '-', 'LineWidth', 1.2, 'Color', [1 1 1]);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'Mixture'; %#ok<AGROW>
    end
    if ~isempty(qGas)
        h = plot(ax, xGas, qGas, '--', 'LineWidth', 1.0, 'Color', [0.9 0.5 0.5]);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'Steam'; %#ok<AGROW>
    end
    if ~isempty(qLiq)
        h = plot(ax, xLiq, qLiq, '--', 'LineWidth', 1.0, 'Color', [0.5 0.8 1.0]);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'Liquid'; %#ok<AGROW>
    end

    if isempty(legendHandles)
        text(ax, 0.5, 0.5, 'No flow data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
    else
        leg = legend(ax, legendHandles, legendLabels, 'Location', 'best');
        applyLegendStyle(leg);
    end

    ax.XLabel.String = 'Position [m]';
    ax.YLabel.String = 'Mass flow [kg/s]';
    ax.Title.String = sprintf('Flow Rates  t = %.3g %s', ...
        profile.timeDays * viewerData.timeScale, viewerData.timeUnitLabel);
    grid(ax, 'on');
    hold(ax, 'off');
end

function renderPrecipitationRate(ax, viewerData, timeIdx, depthIdx)
    cla(ax);
    hold(ax, 'on');

    if isempty(viewerData.scaleThickness) || isempty(viewerData.scaleFractions) || ...
            isempty(viewerData.scaleDensity) || isempty(viewerData.scaleWallArea)
        text(ax, 0.5, 0.5, 'No scale data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
        ax.Title.String = 'Scale Accumulation Rate';
        ax.XLabel.String = sprintf('Time [%s]', viewerData.timeUnitLabel);
        ax.YLabel.String = 'dM/dt [kg/s]';
        hold(ax, 'off');
        return;
    end

    nDepth = size(viewerData.scaleThickness, 1);
    nTimes = size(viewerData.scaleThickness, 2);
    depthIdx = max(1, min(nDepth, round(depthIdx)));
    timeIdx = max(1, min(nTimes, round(timeIdx)));
    [~, rateKgps] = computeScaleMassRateSeries(viewerData, depthIdx);
    if isempty(rateKgps) || ~any(isfinite(rateKgps(2:end)))
        text(ax, 0.5, 0.5, 'No scale data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
        ax.Title.String = sprintf('Scale Accumulation Rate  depth %.1f m', viewerData.x(depthIdx));
        ax.XLabel.String = sprintf('Time [%s]', viewerData.timeUnitLabel);
        ax.YLabel.String = 'dM/dt [kg/s]';
        hold(ax, 'off');
        return;
    end

    tDisplay = viewerData.times(:) * viewerData.timeScale;
    plot(ax, tDisplay, zeros(size(tDisplay)), '--', 'LineWidth', 1.0, ...
        'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    plot(ax, tDisplay, rateKgps, '-', 'LineWidth', 1.2, 'Color', [0.95 0.7 0.15]);

    ax.Title.String = sprintf('Scale Accumulation Rate  depth %.1f m', viewerData.x(depthIdx));
    ax.XLabel.String = sprintf('Time [%s]', viewerData.timeUnitLabel);
    ax.YLabel.String = 'dM/dt [kg/s]';
    grid(ax, 'on');
    hold(ax, 'off');
end

function [massSeriesKg, rateKgps] = computeScaleMassRateSeries(viewerData, depthIdx)
    nTimes = size(viewerData.scaleThickness, 2);
    massSeriesKg = zeros(nTimes, 1);
    rateKgps = NaN(nTimes, 1);
    if nTimes == 0
        return;
    end

    tol = 1e-9;
    depthIdx = max(1, min(size(viewerData.scaleThickness, 1), round(depthIdx)));
    fracDepthIdx = max(1, min(size(viewerData.scaleFractions, 2), round(depthIdx)));
    areaWall = viewerData.scaleWallArea(depthIdx);
    densities = viewerData.scaleDensity(:);
    thicknessSeries = max(viewerData.scaleThickness(depthIdx, :), 0);
    if ~any(thicknessSeries > tol)
        return;
    end
    fractions = reshape(viewerData.scaleFractions(:, fracDepthIdx, :), size(viewerData.scaleFractions, 1), nTimes);

    layerThickness = zeros(1, 0);
    layerFractions = zeros(numel(densities), 0);
    totalThickness = 0;
    arealMass = 0;
    previousFractions = zeros(numel(densities), 1);

    for k = 1:nTimes
        targetThickness = thicknessSeries(k);
        frac = normalizeScaleFractions(fractions(:, k), previousFractions);

        if targetThickness < totalThickness - tol
            amountRemoved = totalThickness - targetThickness;
            while amountRemoved > tol && ~isempty(layerThickness)
                topThickness = layerThickness(end);
                topArealDensity = dot(layerFractions(:, end), densities);
                if topThickness <= amountRemoved + tol
                    arealMass = arealMass - topThickness * topArealDensity;
                    amountRemoved = amountRemoved - topThickness;
                    layerThickness(end) = [];
                    layerFractions(:, end) = [];
                else
                    arealMass = arealMass - amountRemoved * topArealDensity;
                    layerThickness(end) = topThickness - amountRemoved;
                    amountRemoved = 0;
                end
            end
            totalThickness = targetThickness;
        elseif targetThickness > totalThickness + tol
            delta = targetThickness - totalThickness;
            layerThickness(end+1) = delta; %#ok<AGROW>
            layerFractions(:, end+1) = frac; %#ok<AGROW>
            arealMass = arealMass + delta * dot(frac, densities);
            totalThickness = totalThickness + delta;
        elseif ~isempty(layerThickness) && any(frac > tol)
            topThickness = layerThickness(end);
            oldArealDensity = dot(layerFractions(:, end), densities);
            newArealDensity = dot(frac, densities);
            arealMass = arealMass + topThickness * (newArealDensity - oldArealDensity);
            layerFractions(:, end) = frac;
        end

        if ~isempty(layerFractions)
            previousFractions = layerFractions(:, end);
        else
            previousFractions = frac;
        end
        massSeriesKg(k) = max(arealMass, 0) * areaWall;
    end

    dtSec = diff(viewerData.times(:)) * 86400;
    validDt = dtSec > 0;
    deltaMass = diff(massSeriesKg);
    validIdx = find(validDt) + 1;
    rateKgps(validIdx) = deltaMass(validDt) ./ dtSec(validDt);
end

function frac = normalizeScaleFractions(frac, previousFractions)
    tol = 1e-9;
    frac = frac(:);
    frac(~isfinite(frac)) = 0;
    frac = max(frac, 0);
    fracSum = sum(frac);
    if fracSum > tol
        frac = frac / fracSum;
    elseif ~isempty(previousFractions)
        frac = previousFractions(:);
    end
end

function [xMix, qMix] = sumPhaseFlows(profile, x, xGas, qGas, xLiq, qLiq)
    xMix = [];
    qMix = [];

    if ~isempty(profile.Q_mix)
        [xMix, qMix] = alignVectors(x, profile.Q_mix);
        if ~isempty(qMix)
            return;
        end
    end

    if ~isempty(profile.Qm_v) && ~isempty(profile.Qm_l)
        [xV, qV] = alignVectors(x, profile.Qm_v);
        [~, qL] = alignVectors(x, profile.Qm_l);
        mLen = min(numel(qV), numel(qL));
        if mLen > 0
            xMix = xV(1:mLen);
            qMix = qV(1:mLen) + qL(1:mLen);
            return;
        end
    end

    if ~isempty(qGas) && ~isempty(qLiq)
        mLen = min(numel(qGas), numel(qLiq));
        if mLen > 0
            xMix = xGas(1:mLen);
            qMix = qGas(1:mLen) + qLiq(1:mLen);
            return;
        end
    end
end

function tf = hasMeasuredSeries(viewerData, fieldName)
    tf = false;
    data = viewerData.(fieldName);
    tf = isstruct(data) && ~isempty(data.x) && ~isempty(data.values);
end

function plotMeasuredProfile(ax, data, color, tag)
    if nargin < 4 || isempty(tag)
        tag = 'Measured';
    end
    if nargin < 3 || isempty(color)
        color = [0.95 0.8 0.1];
    end

    marker = findobj(ax, 'Type', 'Line', 'Tag', tag);
    if ~isempty(marker)
        set(marker, 'XData', data.x, 'YData', data.value);
        return;
    end
    plot(ax, data.x, data.value, 'o', ...
        'LineStyle', 'none', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'MarkerSize', 6, ...
        'Tag', tag);
end

function [handles, labels] = plotMeasuredSeries(ax, data, colors, tagPrefix)
    handles = [];
    labels = {};
    if nargin < 3 || isempty(colors)
        colors = {};
    end
    if nargin < 4 || isempty(tagPrefix)
        tagPrefix = 'Measured';
    end
    if isempty(data) || isempty(data.values)
        return;
    end

    x = data.x(:).';
    names = data.names;
    if isempty(names)
        names = compose("Series %d", 1:numel(data.values));
    end
    for i = 1:numel(data.values)
        y = data.values{i};
        [xPlot, yPlot] = alignVectors(x, y);
        if isempty(xPlot) || isempty(yPlot)
            continue;
        end
        if ~isempty(colors) && numel(colors) >= i
            color = colors{i};
        else
            color = [0.95 0.8 0.1];
        end
        tag = sprintf('%s_%d', tagPrefix, i);
        h = findobj(ax, 'Type', 'Line', 'Tag', tag);
        if ~isempty(h)
            set(h, 'XData', xPlot, 'YData', yPlot);
        else
            h = plot(ax, xPlot, yPlot, 'o', ...
                'LineStyle', 'none', ...
                'MarkerFaceColor', color, ...
                'MarkerEdgeColor', [0.2 0.2 0.2], ...
                'MarkerSize', 5, ...
                'Tag', tag);
        end
        handles(end+1) = h; %#ok<AGROW>
        labels{end+1} = char(names(min(i, numel(names)))); %#ok<AGROW>
    end
end

function [xOut, dataOut] = alignVectors(x, data)
    if isempty(x) || isempty(data)
        xOut = [];
        dataOut = [];
        return;
    end
    n = min(numel(x), numel(data));
    xOut = x(1:n);
    dataOut = data(1:n);
end

function applyLegendStyle(leg)
    if isempty(leg) || ~isvalid(leg)
        return;
    end
    if isprop(leg, 'ItemTokenSize')
        leg.ItemTokenSize = [12 8];
    end
end

function renderScaleCompositionPlot(ax, viewerData, depthIdx, timeIdx)
    cla(ax);
    hold(ax, 'on');

    state = computeScaleLayerState(viewerData, depthIdx, timeIdx);
    if ~state.hasData
        text(ax, 0.5, 0.5, state.message, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        ax.XDir = 'reverse';
        ax.Title.String = state.titleStr;
        legend(ax, 'off');
        hold(ax, 'off');
        return;
    end

    h = area(ax, state.normThickness, state.areaData);
    palette = getMineralPalette(viewerData);
    if isempty(palette)
        palette = lines(max(numel(h), 1));
    end
    for k = 1:numel(h)
        h(k).FaceColor = palette(min(k, size(palette, 1)), :);
        h(k).EdgeColor = 'none';
    end

    ax.XLim = [0 1];
    ax.YLim = [0 1];
    ax.XDir = 'reverse';
    ax.XLabel.String = 'Normalized thickness [-] (Interior -> Wall)';
    ax.YLabel.String = 'Mineral fraction [-]';
    ax.Title.String = state.titleStr;

    mineralNames = viewerData.scaleMineralNames;
    if ~isempty(mineralNames)
        mineralLabels = cellstr(mineralNames(:));
        leg = legend(ax, mineralLabels, ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal', ...
            'NumColumns', max(1, ceil(numel(mineralLabels) / 2)));
        applyLegendStyle(leg);
    else
        legend(ax, 'off');
    end
    hold(ax, 'off');
end

function renderScaleThicknessHistory(ax, viewerData, depthIdx, timeIdx)
    cla(ax);
    hold(ax, 'on');

    if isempty(viewerData.scaleThickness)
        text(ax, 0.5, 0.5, 'No scale data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        hold(ax, 'off');
        return;
    end

    nDepth = size(viewerData.scaleThickness, 1);
    nTimes = size(viewerData.scaleThickness, 2);
    if nTimes == 0
        text(ax, 0.5, 0.5, 'No scale data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        hold(ax, 'off');
        return;
    end

    depthIdx = max(1, min(nDepth, round(depthIdx)));
    timeIdx = max(1, min(nTimes, round(timeIdx)));
    thicknessSeries = viewerData.scaleThickness(depthIdx, :);
    if ~any(thicknessSeries > 0)
        text(ax, 0.5, 0.5, 'No scale data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.8 0.8 0.8]);
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        ax.Title.String = sprintf('Scale Thickness History  depth %.1f m', viewerData.x(depthIdx));
        hold(ax, 'off');
        return;
    end

    timesDisplay = viewerData.times(:).' * viewerData.timeScale;
    plot(ax, timesDisplay, thicknessSeries, '-', 'LineWidth', 1.2, 'Color', [0.2 0.8 1.0]);
    plot(ax, timesDisplay(timeIdx), thicknessSeries(timeIdx), 'o', ...
        'MarkerFaceColor', [1.0 0.9 0.6], ...
        'MarkerEdgeColor', [0.4 0.3 0.1], ...
        'MarkerSize', 6, ...
        'HandleVisibility', 'off');

    yMax = max(thicknessSeries);
    if yMax <= 0
        yMax = 1;
    end
    ax.YLim = [0 yMax * 1.05];
    ax.XLim = [timesDisplay(1) timesDisplay(end)];
    ax.Title.String = sprintf('Scale Thickness History  depth %.1f m', viewerData.x(depthIdx));
    ax.XLabel.String = sprintf('Time [%s]', viewerData.timeUnitLabel);
    ax.YLabel.String = 'Thickness [m]';
    grid(ax, 'on');
    hold(ax, 'off');
end

function state = computeScaleLayerState(viewerData, depthIdx, timeIdx)
    state = struct('hasData', false, 'titleStr', '', 'message', 'No scale data', ...
        'normThickness', [], 'areaData', [], 'layerThickness', [], ...
        'layerFractions', [], 'totalThickness', 0);

    tol = 1e-9;
    if isempty(viewerData) || isempty(viewerData.x)
        return;
    end

    nDepthPositions = numel(viewerData.x);
    depthInput = round(depthIdx);
    depthForTitle = max(1, min(nDepthPositions, depthInput));
    state.titleStr = sprintf('Scale Composition  depth %.1f m', viewerData.x(depthForTitle));

    if isempty(viewerData.scaleFractions) || isempty(viewerData.scaleThickness)
        return;
    end

    nDepthThickness = size(viewerData.scaleThickness, 1);
    nTimes = size(viewerData.scaleThickness, 2);
    if nDepthThickness == 0 || nTimes == 0
        return;
    end

    depthIdxData = max(1, min(nDepthThickness, depthInput));
    timeIdx = max(1, min(nTimes, round(timeIdx)));
    thicknessSeries = viewerData.scaleThickness(depthIdxData, 1:timeIdx);
    if ~any(thicknessSeries > 0)
        return;
    end

    nMinerals = size(viewerData.scaleFractions, 1);
    nDepthFractions = size(viewerData.scaleFractions, 2);
    if nMinerals == 0
        state.message = 'No mineral data';
        return;
    end

    depthIdxFrac = max(1, min(nDepthFractions, depthInput));
    fractions = reshape(viewerData.scaleFractions(:, depthIdxFrac, 1:timeIdx), nMinerals, timeIdx);
    thicknessSeries = max(thicknessSeries(:).', 0);
    if ~any(thicknessSeries > tol)
        return;
    end

    layerThickness = zeros(1, 0);
    layerFractions = zeros(nMinerals, 0);
    totalThickness = 0;

    for k = 1:numel(thicknessSeries)
        targetThickness = thicknessSeries(k);
        frac = fractions(:, k);
        frac(~isfinite(frac)) = 0;
        frac = max(frac, 0);
        fracSum = sum(frac);
        if fracSum > tol
            frac = frac / fracSum;
        elseif ~isempty(layerFractions)
            frac = layerFractions(:, end);
        else
            frac(:) = 0;
        end

        if targetThickness < totalThickness - tol
            amountRemoved = totalThickness - targetThickness;
            while amountRemoved > tol && ~isempty(layerThickness)
                topThickness = layerThickness(end);
                if topThickness <= amountRemoved + tol
                    amountRemoved = amountRemoved - topThickness;
                    layerThickness(end) = [];
                    layerFractions(:, end) = [];
                else
                    layerThickness(end) = topThickness - amountRemoved;
                    amountRemoved = 0;
                end
            end
            totalThickness = targetThickness;
        elseif targetThickness > totalThickness + tol
            delta = targetThickness - totalThickness;
            layerThickness(end+1) = delta; %#ok<AGROW>
            layerFractions(:, end+1) = frac; %#ok<AGROW>
            totalThickness = totalThickness + delta;
        else
            if ~isempty(layerFractions) && any(frac > tol)
                layerFractions(:, end) = frac;
            end
            totalThickness = targetThickness;
        end
    end

    if isempty(layerThickness)
        return;
    end

    totalThickness = sum(layerThickness);
    if totalThickness <= tol
        return;
    end

    nLayers = numel(layerThickness);
    xVals = zeros(2 * nLayers, 1);
    yVals = zeros(2 * nLayers, nMinerals);
    cursor = 0;
    idx = 1;
    for k = 1:nLayers
        segThickness = layerThickness(k);
        if segThickness <= tol
            continue;
        end
        frac = layerFractions(:, k).';
        xVals(idx) = cursor;
        yVals(idx, :) = frac;
        idx = idx + 1;
        cursor = cursor + segThickness;
        xVals(idx) = cursor;
        yVals(idx, :) = frac;
        idx = idx + 1;
    end

    xVals = xVals(1:idx-1);
    yVals = yVals(1:idx-1, :);
    if numel(xVals) < 2 || cursor <= tol
        return;
    end

    normThickness = xVals / cursor;
    normThickness(end) = 1;
    totalFrac = sum(yVals, 2);
    nz = totalFrac > tol;
    yVals(nz, :) = yVals(nz, :) ./ totalFrac(nz);
    yVals(~nz, :) = 0;

    state.hasData = true;
    state.message = '';
    state.normThickness = normThickness;
    state.areaData = yVals;
    state.layerThickness = layerThickness;
    state.layerFractions = layerFractions;
    state.totalThickness = cursor;
end

function fracSlice = computeScaleFractionSlice(viewerData, timeIdx)
    if isempty(viewerData) || isempty(viewerData.scaleFractions)
        fracSlice = [];
        return;
    end

    nMinerals = size(viewerData.scaleFractions, 1);
    nDepth = numel(viewerData.x);
    fracSlice = zeros(nMinerals, nDepth);
    tol = 1e-9;
    for depthIdx = 1:nDepth
        state = computeScaleLayerState(viewerData, depthIdx, timeIdx);
        if ~state.hasData || isempty(state.layerThickness)
            continue;
        end
        contributions = state.layerFractions .* state.layerThickness;
        mineralTotals = sum(contributions, 2);
        total = sum(mineralTotals);
        if total > tol
            fracSlice(:, depthIdx) = mineralTotals / total;
        end
    end
end

function plotGeometryScaleFill(ax, depthPositions, baseDp, currDp, fracSlice, palette)
    tol = 1e-9;
    defaultColor = [0.6 0.6 0.6];
    n = numel(depthPositions);
    if n < 2 || isempty(baseDp) || isempty(currDp)
        return;
    end

    baseDp = alignLengthLocal(baseDp, n);
    currDp = alignLengthLocal(currDp, n);
    if nargin < 5 || isempty(fracSlice)
        fracSlice = [];
    else
        if size(fracSlice, 2) < n
            fracSlice(:, end+1:n) = 0;
        elseif size(fracSlice, 2) > n
            fracSlice = fracSlice(:, 1:n);
        end
    end

    nodeColors = repmat(defaultColor, n, 1);
    nodeHasData = false(n, 1);
    if ~isempty(fracSlice)
        for i = 1:n
            nodeColors(i, :) = mineralMixColor(fracSlice(:, i), palette, defaultColor);
            nodeHasData(i) = any(fracSlice(:, i) > tol);
        end
    end

    for i = 1:n-1
        deltaTop = baseDp(i) - currDp(i);
        deltaBot = baseDp(i+1) - currDp(i+1);
        if deltaTop <= tol && deltaBot <= tol
            continue;
        end
        outerTop = max(baseDp(i), currDp(i));
        innerTop = min(baseDp(i), currDp(i));
        outerBot = max(baseDp(i+1), currDp(i+1));
        innerBot = min(baseDp(i+1), currDp(i+1));
        if outerTop - innerTop <= tol && outerBot - innerBot <= tol
            continue;
        end

        colorTop = nodeColors(i, :);
        colorBot = nodeColors(i+1, :);
        if nodeHasData(i) && ~nodeHasData(i+1)
            segColor = colorTop;
        elseif nodeHasData(i+1) && ~nodeHasData(i)
            segColor = colorBot;
        else
            segColor = 0.5 * (colorTop + colorBot);
        end
        segColor = min(max(segColor, 0), 1);

        patch(ax, ...
            [innerTop, outerTop, outerBot, innerBot], ...
            [depthPositions(i), depthPositions(i), depthPositions(i+1), depthPositions(i+1)], ...
            segColor, 'EdgeColor', 'none', 'FaceAlpha', 1.0, 'Tag', 'ScaleDepositFill');
    end
end

function color = mineralMixColor(fracVec, palette, defaultColor)
    if nargin < 3 || isempty(defaultColor)
        defaultColor = [0.6 0.6 0.6];
    end
    if isempty(fracVec)
        color = defaultColor;
        return;
    end
    weights = max(fracVec(:), 0);
    total = sum(weights);
    if total <= 1e-12 || isempty(palette)
        color = defaultColor;
        return;
    end
    weights = weights / total;
    if size(palette, 1) < numel(weights)
        palette = [palette; repmat(palette(end, :), numel(weights) - size(palette, 1), 1)];
    end
    color = weights.' * palette(1:numel(weights), :);
    color = min(max(color, 0), 1);
end

function palette = getMineralPalette(viewerData)
    nMinerals = 0;
    if ~isempty(viewerData.scaleFractions)
        nMinerals = size(viewerData.scaleFractions, 1);
    elseif ~isempty(viewerData.scaleMineralNames)
        nMinerals = numel(viewerData.scaleMineralNames);
    end
    if nMinerals > 0
        palette = lines(nMinerals);
    else
        palette = [];
    end
end

function out = alignLengthLocal(data, n)
    data = data(:);
    if numel(data) >= n
        out = data(1:n);
    elseif isempty(data)
        out = zeros(n, 1);
    else
        out = [data; repmat(data(end), n - numel(data), 1)];
    end
end
