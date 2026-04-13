function plotResultsOnAxes(axesStruct, state, Y, titleText)
% Plot solution Y onto existing axes in axesStruct, keeping previous lines.
% This is for transient overlays (no clearing).

x = state.x;

% Extract solution variables
pressure = Y(1, :);
enthalpy = Y(2, :);
velocity = Y(3, :);

% Derived quantities
n = length(x);
[alpha_g, alpha_l, rho_g, rho_l, ~, ~, ~, temperature] = ...
    calculatePhaseProperties(pressure, enthalpy, state);
thetaNodes = state.gravityThetaNode(:).';
[u_g, u_l] = calculatePhaseVelocities(velocity, alpha_g, rho_g, rho_l, state.Dp(:).', temperature, thetaNodes, state);

% Mineral palette shared between panels
mineralNames = string(state.chem.mineralNames(:));
nMinerals = numel(mineralNames);
mineralColors = zeros(0, 3);
if nMinerals > 0
    mineralColors = lines(nMinerals);
end

% Colors for overlay lines
c_pressure = [0.0 1.0 1.0];
c_temp     = [1.0 0.0 0.0];
c_u_mix    = [0.2 0.8 0.2];
c_u_g      = [1.0 0.3 0.3];
c_u_l      = [0.3 0.9 0.9];
c_q_mix    = [1.0 1.0 1.0];
c_q_g      = [0.9 0.5 0.5];
c_q_l      = [0.5 0.8 1.0];
c_feed_p   = [0.95 0.9 0.2];
c_feed_t   = [1.0 0.6 0.2];

% Pressure
pScale = state.pressureUnitScale;
pLabel = char(state.pressureUnitLabel);
[feedXP, feedP] = getFeedzonePressurePlot(state, pScale);
[feedXT, feedT] = getFeedzoneTemperaturePlot(state);
if ~isempty(state.measured.pressure.position) && ~isempty(state.measured.pressure.value)
    xMeas = getMeasuredPositions(state, state.measured.pressure);
    plotOrUpdateMeasured(axesStruct.ax1, xMeas, state.measured.pressure.value, ...
        'MeasuredPressure', [0.95 0.8 0.1], [0.95 0.6 0.0]);
end
initP = state.initProfiles.pressure;
if ~isempty(initP.position) && ~isempty(initP.value)
    xInit = getMeasuredPositions(state, initP);
    [xInit, pInit] = alignSeries(xInit, initP.value / pScale);
    plotOrUpdateProfileLine(axesStruct.ax1, xInit, pInit, ...
        'InitPressure', '--', [0.6 0.9 0.9], 1.0);
end
plot(axesStruct.ax1, x, pressure / pScale, '-', 'LineWidth', 1.0, 'Color', c_pressure, ...
    'DisplayName', 'Pressure');
if ~isempty(feedP)
    plotOrUpdateMarkers(axesStruct.ax1, feedXP, feedP, 'FeedzonePressure', '*', c_feed_p, 5);
end
try
    ylabel(axesStruct.ax1, sprintf('P [%s]', pLabel));
end
if ~isempty(titleText), axesStruct.ax1.Title.String = ['Pressure  ' titleText]; end

% Temperature (C)
if ~isempty(state.measured.temperature.position) && ~isempty(state.measured.temperature.value)
    xMeas = getMeasuredPositions(state, state.measured.temperature);
    plotOrUpdateMeasured(axesStruct.ax2, xMeas, state.measured.temperature.value, ...
        'MeasuredTemperature', [0.4 0.8 1.0], [0.1 0.5 0.9]);
end
tPlotMin = temperature(end) - 273.15 - 10;
tPlotMax = max(temperature) - 273.15 + 10;
initT = state.initProfiles.temperature;
if ~isempty(initT.position) && ~isempty(initT.value)
    xInit = getMeasuredPositions(state, initT);
    [xInit, tInit] = alignSeries(xInit, initT.value);
    tInit = min(max(tInit, tPlotMin), tPlotMax);
    plotOrUpdateProfileLine(axesStruct.ax2, xInit, tInit, ...
        'InitTemperature', '--', [1.0 0.5 0.5], 1.0);
end
if ~isempty(state.measured.pressure.position) && ~isempty(state.measured.pressure.value)
    xSat = getMeasuredPositions(state, state.measured.pressure);
    [xSat, pMeas] = alignSeries(xSat, state.measured.pressure.value);
    if ~isempty(pMeas)
        tSatC = [];
        try
            pPa = pMeas * pScale;
            tSatC = IAPWS_IF97('Tsat_p', pPa / 1e6) - 273.15;
        catch
            tSatC = [];
        end
        if ~isempty(tSatC)
            tSatC = min(max(tSatC, tPlotMin), tPlotMax);
            plotOrUpdateProfileLine(axesStruct.ax2, xSat, tSatC, ...
                'TsatFromPressure', ':', [0.7 0.7 0.9], 1.0);
        end
    end
end
initP = state.initProfiles.pressure;
if ~isempty(initP.position) && ~isempty(initP.value)
    xSatInit = getMeasuredPositions(state, initP);
    [xSatInit, pInitPa] = alignSeries(xSatInit, initP.value);
    if ~isempty(pInitPa)
        tSatInitC = [];
        try
            tSatInitC = IAPWS_IF97('Tsat_p', pInitPa / 1e6) - 273.15;
        catch
            tSatInitC = [];
        end
        if ~isempty(tSatInitC)
            tSatInitC = min(max(tSatInitC, tPlotMin), tPlotMax);
            plotOrUpdateProfileLine(axesStruct.ax2, xSatInit, tSatInitC, ...
                'TsatFromInitPressure', ':', [0.6 0.6 0.6], 1.0);
        end
    end
end
plot(axesStruct.ax2, x, temperature-273.15, '-', 'LineWidth', 1.0, 'Color', c_temp, ...
    'DisplayName', 'Temperature');
if ~isempty(feedT)
    plotOrUpdateMarkers(axesStruct.ax2, feedXT, feedT, 'FeedzoneTemperature', '*', c_feed_t, 5);
end
if ~isempty(titleText), axesStruct.ax2.Title.String = ['Temperature  ' titleText]; end
try
    xlim(axesStruct.ax2, [0, state.x(end)]);
catch
end

% Velocities
legendHandles = [];
legendLabels = {};
mh = [];
mlabels = {};
if ~isempty(state.measured.velocity.values)
    measColors = {[0.95 0.7 0.1], [0.65 0.6 1.0], [0.9 0.9 0.9]};
    [mh, mlabels] = plotMeasuredMultiSeries(axesStruct.ax3, state, state.measured.velocity, ...
        'MeasuredVel', measColors, []);
end
legendHandles(end+1) = plot(axesStruct.ax3, x, velocity, '-', 'LineWidth', 1.0, ...
    'Color', c_u_mix, 'DisplayName', 'Mixture');
legendLabels{end+1} = 'Mixture';
legendHandles(end+1) = plot(axesStruct.ax3, x, u_g, '--', 'LineWidth', 0.8, ...
    'Color', c_u_g, 'DisplayName', 'Gas');
legendLabels{end+1} = 'Gas';
legendHandles(end+1) = plot(axesStruct.ax3, x, u_l, '--', 'LineWidth', 0.8, ...
    'Color', c_u_l, 'DisplayName', 'Liquid');
legendLabels{end+1} = 'Liquid';
if ~isempty(mh)
    legendHandles = [legendHandles, mh];
    legendLabels = [legendLabels, mlabels];
end
if ~isempty(legendHandles)
    leg = legend(axesStruct.ax3, legendHandles, legendLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
    applyLegendTokenSize(leg);
else
    legend(axesStruct.ax3, 'off');
end
if ~isempty(titleText), axesStruct.ax3.Title.String = ['Velocities  ' titleText]; end

% Mineral saturation indices (bottom-left)
try
    axSI = axesStruct.ax4;
    cla(axSI);
    hold(axSI, 'on');
    plot(axSI, x, zeros(size(x)), '--', 'LineWidth', 1.0, 'Color', [0.6 0.6 0.6]);
    siMatrix = state.chem.saturationIndices;
    siNames = string(state.chem.mineralNames(:));
    nMin = min(size(siMatrix, 1), numel(siNames));
    cols = min(size(siMatrix, 2), numel(x));
    colors = mineralColors;
    if size(colors, 1) < nMin
        colors = lines(max(nMin, 1));
    end
    colors = colors(1:nMin, :);
    for idx = 1:nMin
        plot(axSI, x(1:cols), siMatrix(idx, 1:cols), '-', 'LineWidth', 1.0, ...
            'Color', colors(idx, :), 'DisplayName', char(siNames(idx)));
    end
    hold(axSI, 'off');
    if nMin > 0
        leg = legend(axSI, cellstr(siNames(1:nMin)), 'Location', 'southoutside', ...
            'Orientation', 'horizontal', 'NumColumns', max(1, ceil(nMin/2)));
        applyLegendTokenSize(leg);
    else
        legend(axSI, 'off');
    end
    ylabel(axSI, 'SI [-]');
    ylim(axSI, [-2 2]);
    if ~isempty(titleText)
        axSI.Title.String = ['Saturation Index  ' titleText];
    end
catch ME
    warning('Saturation index plot failed: %s', ME.message);
end

% Element mass concentration plot (bottom-middle)
axElem = axesStruct.ax5;
cla(axElem);
axElem.YScale = 'log';
axElem.XMinorGrid = 'off';
axElem.YMinorGrid = 'off';
hold(axElem, 'on');
speciesIdx = state.chem.elementSpeciesIndex(:);
plotLabels = string(state.chem.plotLabels(:));
gasMask = state.chem.gasMask(:);
nElems = numel(speciesIdx);
cols = min(size(state.C, 2), numel(x));
chemColors = chemPlotPalette(numel(plotLabels));
chemPpm = liquidMassfracToPpmV2(state.C(:, 1:cols));
hElem = gobjects(nElems, 1);
for elemIdx = 1:nElems
    spIdx = speciesIdx(elemIdx);
    lineStyle = '-';
    if gasMask(spIdx)
        lineStyle = '--';
    end
    hElem(elemIdx) = semilogy(axElem, x(1:cols), chemPpm(spIdx, 1:cols), ...
        'LineStyle', lineStyle, 'LineWidth', 1.2, 'Color', chemColors(spIdx, :), ...
        'DisplayName', char(plotLabels(spIdx)));
end
hold(axElem, 'off');

if nElems > 0
    leg = legend(axElem, hElem, cellstr(plotLabels(speciesIdx)), 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'NumColumns', max(1, ceil(nElems/3)));
    applyLegendTokenSize(leg);
else
    legend(axElem, 'off');
end
ylim(axElem, [1e-2 max(1e-1, max(chemPpm(:)))]);
ylabel(axElem, 'c [ppm]');
if ~isempty(titleText)
    axElem.Title.String = ['Mass Concentration  ' titleText];
end

% Flow rates (bottom-right of left block)
legendHandles = [];
legendLabels = {};
stateFlux = refreshHydroFluxCache(state, Y);
qGasFull = stateFlux.Q_v(1:n);
qLiqFull = stateFlux.Q_l(1:n);
[xMix, qMix] = alignSeries(x, qGasFull + qLiqFull);
[xGas, qGas] = alignSeries(x, qGasFull);
[xLiq, qLiq] = alignSeries(x, qLiqFull);
if exist('qMix', 'var') && ~isempty(qMix)
    legendHandles(end+1) = plot(axesStruct.ax6, xMix, qMix, '-', 'LineWidth', 1.0, ...
        'Color', c_q_mix, 'DisplayName', 'Total');
    legendLabels{end+1} = 'Total';
end
if exist('qGas', 'var') && ~isempty(qGas)
    legendHandles(end+1) = plot(axesStruct.ax6, xGas, qGas, '--', 'LineWidth', 0.8, ...
        'Color', c_q_g, 'DisplayName', 'Steam');
    legendLabels{end+1} = 'Steam';
end
if exist('qLiq', 'var') && ~isempty(qLiq)
    legendHandles(end+1) = plot(axesStruct.ax6, xLiq, qLiq, '--', 'LineWidth', 0.8, ...
        'Color', c_q_l, 'DisplayName', 'Liquid');
    legendLabels{end+1} = 'Liquid';
end
if ~isempty(legendHandles)
    leg = legend(axesStruct.ax6, legendHandles, legendLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
    applyLegendTokenSize(leg);
else
    legend(axesStruct.ax6, 'off');
end
if ~isempty(titleText), axesStruct.ax6.Title.String = ['Flow Rates  ' titleText]; end
try
    ylabel(axesStruct.ax6, 'Mass flow [kg/s]');
end

%% Geometry column at right: plot diameters and color scale deposits by composition
try
    axGeom = axesStruct.axGeom;
    xx = state.x(:);
    Dp  = state.Dp(:);
    Dp0 = state.Dp0(:);
    scaleDp = [];
    scaleX = [];
    if isfield(state, 'scaleProfile') && ~isempty(state.scaleProfile) && ...
            isfield(state.scaleProfile, 'position') && isfield(state.scaleProfile, 'diameter')
        scaleX = state.scaleProfile.position(:);
        scaleDp = state.scaleProfile.diameter(:);
        [scaleX, scaleDp] = alignSeries(scaleX, scaleDp);
    end
    geomCount = min([numel(xx), numel(Dp), numel(Dp0)]);
    xx = xx(1:geomCount);
    Dp = Dp(1:geomCount);
    Dp0 = Dp0(1:geomCount);

    if geomCount == 0
        if ~isempty(titleText), axGeom.Title.String = ['Geometry  ' titleText]; end
    else
        % Prepare mineral thickness and colors per depth slice
        thickness = [];
        currFrac = [];
        palette = [];
        if ~isempty(state.chem.mineralThickness) && size(state.chem.mineralThickness, 1) > 0
            fullThickness = state.chem.mineralThickness;
            nRows = size(fullThickness, 1);
            colorRows = min(3, nRows);
            if colorRows > 0
                thickness = fullThickness(1:colorRows, :);
                if size(thickness, 2) < geomCount
                    deficit = geomCount - size(thickness, 2);
                    thickness = [thickness, zeros(colorRows, deficit)];
                elseif size(thickness, 2) > geomCount
                    thickness = thickness(:, 1:geomCount);
                end
                thickness = max(thickness, 0);

                palette = mineralColors;
                if size(palette, 1) < colorRows
                    palette = lines(max(colorRows, 1));
                end
                palette = palette(1:colorRows, :);

                if ~isempty(state.chem.totalScaleThickness)
                    totalThickness = state.chem.totalScaleThickness(:);
                else
                    totalThickness = sum(fullThickness, 1).';
                end
                if numel(totalThickness) < geomCount
                    totalThickness = [totalThickness; zeros(geomCount - numel(totalThickness), 1)];
                elseif numel(totalThickness) > geomCount
                    totalThickness = totalThickness(1:geomCount);
                end
                denom = max(totalThickness(:).', 1e-12);
                currFrac = bsxfun(@rdivide, thickness, denom);
                zeroMask = totalThickness(:).' <= 1e-12;
                if any(zeroMask)
                    currFrac(:, zeroMask) = 0;
                end
            end
        end

        % Retrieve previous diameter profile for incremental fill
        axData = axGeom.UserData;
        lastDp = axData.lastDp(:);
        prevFrac = axData.lastFractions;
        hasPrev = ~isempty(lastDp);
        if hasPrev
            if numel(lastDp) < geomCount
                lastDp = [lastDp; repmat(lastDp(end), geomCount - numel(lastDp), 1)];
            elseif numel(lastDp) > geomCount
                lastDp = lastDp(1:geomCount);
            end
            if ~isempty(prevFrac)
                if size(prevFrac, 2) < geomCount
                    prevFrac = [prevFrac, zeros(size(prevFrac,1), geomCount - size(prevFrac,2))];
                elseif size(prevFrac, 2) > geomCount
                    prevFrac = prevFrac(:, 1:geomCount);
                end
            end
        end

        % Remove existing scale segments and redraw per depth slice
        delete(findall(axGeom, 'Tag', 'ScaleSegment'));
        tol = 1e-9;
        if hasPrev && geomCount > 1
            if isempty(currFrac)
                currFrac = zeros(0, geomCount);
            end
            if isempty(prevFrac)
                prevFrac = zeros(size(currFrac));
            elseif size(prevFrac,1) ~= size(currFrac,1)
                prevFrac = zeros(size(currFrac));
            end
            for i = 1:geomCount-1
                dPrevTop = lastDp(i);
                dPrevBot = lastDp(i+1);
                dCurrTop = Dp(i);
                dCurrBot = Dp(i+1);
                % Old scale portion (previous step) between pristine and last diameter
                oldOuterTop = max(Dp0(i), dPrevTop);
                oldInnerTop = min(Dp0(i), dPrevTop);
                oldOuterBot = max(Dp0(i+1), dPrevBot);
                oldInnerBot = min(Dp0(i+1), dPrevBot);
                if (oldOuterTop - oldInnerTop > tol) || (oldOuterBot - oldInnerBot > tol)
                    oldColor = [0.6 0.6 0.6];
                    if ~isempty(prevFrac) && ~isempty(palette)
                        prevTop = prevFrac(:, i);
                        prevBot = prevFrac(:, i+1);
                        if sum(prevTop) > 1e-12
                            colorTopPrev = min(max(prevTop.' * palette, 0), 1);
                        else
                            colorTopPrev = oldColor;
                        end
                        if sum(prevBot) > 1e-12
                            colorBotPrev = min(max(prevBot.' * palette, 0), 1);
                        else
                            colorBotPrev = oldColor;
                        end
                        oldColor = 0.5 * (colorTopPrev + colorBotPrev);
                    end
                    xOld = [oldInnerTop, oldOuterTop, oldOuterBot, oldInnerBot];
                    yOld = [xx(i),        xx(i),      xx(i+1),   xx(i+1)];
                    patch(axGeom, xOld, yOld, oldColor, 'EdgeColor', 'none', ...
                        'FaceAlpha', 0.55, 'Tag', 'ScaleSegment');
                end

                if abs(dPrevTop - dCurrTop) < tol && abs(dPrevBot - dCurrBot) < tol
                    continue;
                end

                if (dCurrTop >= dPrevTop - tol) && (dCurrBot >= dPrevBot - tol)
                    continue; % no new scale added in this slice
                end

                outerTop = max(dPrevTop, dCurrTop);
                innerTop = min(dPrevTop, dCurrTop);
                outerBot = max(dPrevBot, dCurrBot);
                innerBot = min(dPrevBot, dCurrBot);
                if outerTop - innerTop < tol && outerBot - innerBot < tol
                    continue;
                end

                % Newly deposited scale between last and current diameter
                segColor = [0.6 0.6 0.6];
                if ~isempty(currFrac) && ~isempty(palette)
                    deltaTop = max(currFrac(:, i) - prevFrac(:, i), 0);
                    deltaBot = max(currFrac(:, i+1) - prevFrac(:, i+1), 0);
                    if sum(deltaTop) > 1e-12
                        deltaTop = deltaTop / sum(deltaTop);
                    elseif sum(currFrac(:, i)) > 1e-12
                        deltaTop = currFrac(:, i) / sum(currFrac(:, i));
                    else
                        deltaTop = [];
                    end
                    if sum(deltaBot) > 1e-12
                        deltaBot = deltaBot / sum(deltaBot);
                    elseif sum(currFrac(:, i+1)) > 1e-12
                        deltaBot = currFrac(:, i+1) / sum(currFrac(:, i+1));
                    else
                        deltaBot = [];
                    end

                    colorTop = segColor;
                    if ~isempty(deltaTop)
                        colorTop = min(max(deltaTop.' * palette, 0), 1);
                    end
                    colorBot = segColor;
                    if ~isempty(deltaBot)
                        colorBot = min(max(deltaBot.' * palette, 0), 1);
                    end

                    segColor = 0.5 * (colorTop + colorBot);
                    segColor = min(max(segColor, 0), 1);
                end

                xPoly = [innerTop, outerTop, outerBot, innerBot];
                yPoly = [xx(i),    xx(i),   xx(i+1),  xx(i+1)];
                patch(axGeom, xPoly, yPoly, segColor, 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.65, 'Tag', 'ScaleSegment');
            end
        end

        % Update diameter curves without duplicating line objects
        currLine = findobj(axGeom, 'Type', 'Line', 'Tag', 'CurrentDiameterLine');
        if numel(currLine) > 1
            delete(currLine(2:end));
            currLine = currLine(1);
        end
        if isempty(currLine)
            currLine = plot(axGeom, Dp, xx, '-', 'LineWidth', 1.4, 'Color', [0.2 0.8 0.2], ...
                'Tag', 'CurrentDiameterLine', 'DisplayName', 'Current diameter');
        else
            set(currLine, 'XData', Dp, 'YData', xx);
        end

        baseLine = findobj(axGeom, 'Type', 'Line', 'Tag', 'PristineDiameterLine');
        if numel(baseLine) > 1
            delete(baseLine(2:end));
            baseLine = baseLine(1);
        end
        if isempty(baseLine)
            baseLine = plot(axGeom, Dp0, xx, '--', 'LineWidth', 1.0, 'Color', [0.7 0.7 0.7], ...
                'Tag', 'PristineDiameterLine', 'DisplayName', 'Pristine diameter');
        else
            set(baseLine, 'XData', Dp0, 'YData', xx);
        end

        scaleLine = findobj(axGeom, 'Type', 'Line', 'Tag', 'ScaledDiameterLine');
        if numel(scaleLine) > 1
            delete(scaleLine(2:end));
            scaleLine = scaleLine(1);
        end
        if ~isempty(scaleDp)
            if isempty(scaleLine)
                scaleLine = plot(axGeom, scaleDp, scaleX, ':', 'LineWidth', 1.1, 'Color', [1.0 0.8 0.2], ...
                    'Tag', 'ScaledDiameterLine', 'DisplayName', 'Scaled diameter');
            else
                set(scaleLine, 'XData', scaleDp, 'YData', scaleX);
            end
        elseif ~isempty(scaleLine)
            delete(scaleLine);
            scaleLine = gobjects(0);
        end

        try
            uistack(baseLine, 'top');
            if ~isempty(scaleLine)
                uistack(scaleLine, 'top');
            end
            uistack(currLine, 'top');
        catch
        end

        allDiameters = [Dp; Dp0];
        if hasPrev
            allDiameters = [allDiameters; lastDp];
        end
        if ~isempty(scaleDp)
            allDiameters = [allDiameters; scaleDp];
        end
        xmin = min(allDiameters);
        xmax = max(allDiameters);
        pad = 0.02 * max(eps, xmax);
        xlim(axGeom, [max(0, xmin - pad), xmax + pad]);
        if ~isempty(titleText)
            axGeom.Title.String = ['Geometry  ' titleText];
        end

        % Store current diameter for next overlay
        axData.lastDp = Dp;
        if ~isempty(currFrac)
            axData.lastFractions = currFrac;
        else
            axData.lastFractions = [];
        end
        axGeom.UserData = axData;
    end
catch ME
    warning('Geometry column plot failed: %s', ME.message);
end
try
    ylim(axesStruct.axGeom, [0, state.x(end)]);
catch
end
end

function plotOrUpdateMeasured(ax, xdata, ydata, tag, edgeColor, faceColor)
% Create or update measured profile markers on an axis.
if isempty(xdata) || isempty(ydata)
    return;
end

xdata = xdata(:)';
ydata = ydata(:)';

h = findobj(ax, 'Type', 'Line', 'Tag', tag);
if isempty(h) || ~isvalid(h)
    plot(ax, xdata, ydata, '.', 'LineStyle', 'none', 'MarkerSize', 6, ...
        'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', faceColor, ...
        'Tag', tag, 'DisplayName', prettifyTag(tag));
else
    set(h, 'XData', xdata, 'YData', ydata);
    if isprop(h, 'DisplayName') && isempty(h.DisplayName)
        h.DisplayName = prettifyTag(tag);
    end
end
end

function plotOrUpdateProfileLine(ax, xdata, ydata, tag, lineStyle, color, lineWidth)
% Create or update dashed profile line overlays.
if isempty(xdata) || isempty(ydata)
    return;
end

xdata = xdata(:)';
ydata = ydata(:)';

h = findobj(ax, 'Type', 'Line', 'Tag', tag);
if isempty(h) || ~isvalid(h)
    plot(ax, xdata, ydata, 'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
        'Color', color, 'Tag', tag, 'DisplayName', prettifyTag(tag));
else
    set(h, 'XData', xdata, 'YData', ydata, 'LineStyle', lineStyle, ...
        'LineWidth', lineWidth, 'Color', color);
    if isprop(h, 'DisplayName') && isempty(h.DisplayName)
        h.DisplayName = prettifyTag(tag);
    end
end
end

function applyLegendTokenSize(leg)
% Shorten legend line tokens for compact layouts.
if isempty(leg) || ~isvalid(leg)
    return;
end
if isprop(leg, 'ItemTokenSize')
    leg.ItemTokenSize = [8 8];
end
end

function plotOrUpdateMarkers(ax, xdata, ydata, tag, marker, color, markerSize)
% Create or update marker-only series on an axis.
if isempty(xdata) || isempty(ydata)
    return;
end

xdata = xdata(:)';
ydata = ydata(:)';

h = findobj(ax, 'Type', 'Line', 'Tag', tag);
if isempty(h) || ~isvalid(h)
    plot(ax, xdata, ydata, 'LineStyle', 'none', 'Marker', marker, ...
        'MarkerSize', markerSize, 'MarkerEdgeColor', color, ...
        'MarkerFaceColor', color, 'Tag', tag, 'DisplayName', prettifyTag(tag));
else
    set(h, 'XData', xdata, 'YData', ydata, 'LineStyle', 'none', 'Marker', marker, ...
        'MarkerSize', markerSize, 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    if isprop(h, 'DisplayName') && isempty(h.DisplayName)
        h.DisplayName = prettifyTag(tag);
    end
end
end

function [xPlot, pPlot] = getFeedzonePressurePlot(state, pScale)
% Return feedzone pressure markers scaled to current pressure units.
xPlot = [];
pPlot = [];
% Prefer plotting all feedzone cells when available
idx = find(state.feedzone_cells(:) == 1);
if ~isempty(idx)
    xBase = state.x(idx);
    pRes = state.feedzone_P_res(idx);
    [xPlot, pPlot] = alignSeries(xBase, pRes / pScale);
    return;
end

fz = state.feedzones;
depth = fz.depth(:);
pRes = fz.P_res(:);
if isempty(depth) || isempty(pRes)
    return;
end
n = min(numel(depth), numel(pRes));
depth = depth(1:n);
pRes = pRes(1:n);
if isempty(depth)
    return;
end
xBase = getFeedzonePositions(state, fz, depth);
[xPlot, pPlot] = alignSeries(xBase, pRes / pScale);
end

function [xPlot, tPlot] = getFeedzoneTemperaturePlot(state)
% Return feedzone temperature markers in Celsius.
xPlot = [];
tPlot = [];
% Prefer plotting all feedzone cells when available
idx = find(state.feedzone_cells(:) == 1);
if ~isempty(idx)
    xBase = state.x(idx);
    pRes = state.feedzone_P_res(idx);
    hRes = state.feedzone_H_res(idx);
    if isempty(xBase)
        return;
    end
    tKelvin = [];
    if state.interp
        try
            tKelvin = state.Temp(hRes, pRes);
        catch
            tKelvin = [];
        end
    end
    if isempty(tKelvin)
        try
            tKelvin = IAPWS_IF97('T_ph', pRes / 1e6, hRes / 1e3);
        catch
            tKelvin = nan(size(pRes));
        end
    end
    [xPlot, tPlot] = alignSeries(xBase, tKelvin - 273.15);
    return;
end

fz = state.feedzones;
depth = fz.depth(:);
pRes = fz.P_res(:);
hRes = fz.H_res(:);
if isempty(depth) || isempty(pRes) || isempty(hRes)
    return;
end
n = min([numel(depth), numel(pRes), numel(hRes)]);
depth = depth(1:n);
pRes = pRes(1:n);
hRes = hRes(1:n);
if isempty(depth)
    return;
end
xBase = getFeedzonePositions(state, fz, depth);
tKelvin = [];
if state.interp
    try
        tKelvin = state.Temp(hRes, pRes);
    catch
        tKelvin = [];
    end
end
if isempty(tKelvin)
    try
        tKelvin = IAPWS_IF97('T_ph', pRes / 1e6, hRes / 1e3);
    catch
        tKelvin = nan(size(pRes));
    end
end
[xPlot, tPlot] = alignSeries(xBase, tKelvin - 273.15);
end

function xBase = getFeedzonePositions(state, fz, depth)
% Convert feedzone depth to simulation position (distance from bottom).
xBase = [];
xBase = fz.position(:);
if isempty(xBase)
    xBase = depth(:);
    xBase = state.x(end) - depth(:);
    xBase = max(0, min(state.Lp, xBase));
end
xBase = xBase(:);
end

function [handles, labels] = plotMeasuredMultiSeries(ax, state, measStruct, tagPrefix, colors, nameFilter)
% Overlay multiple measured series (depth-based) on an axis.
handles = [];
labels = {};

if isempty(measStruct) || isempty(measStruct.values)
    return;
end

xBase = getMeasuredPositions(state, measStruct);
names = measStruct.names;
vals = measStruct.values;
if isempty(names)
    names = strings(1, numel(vals));
    for i = 1:numel(vals)
        names(i) = sprintf('Series %d', i);
    end
end

idxKeep = true(1, numel(vals));
if ~isempty(nameFilter)
    for i = 1:numel(vals)
        try
            idxKeep(i) = nameFilter(names(i));
        catch
            idxKeep(i) = true;
        end
    end
end

for i = 1:numel(vals)
    if ~idxKeep(i)
        continue;
    end
    y = vals{i};
    [xPlot, yPlot] = alignSeries(xBase, y);
    if isempty(xPlot) || isempty(yPlot)
        continue;
    end
    color = [];
    if ~isempty(colors) && numel(colors) >= i
        color = colors{i};
    end
    if isempty(color)
        color = [0.95 0.8 0.1];
    end
    tag = sprintf('%s_%d', tagPrefix, i);
    h = findobj(ax, 'Type', 'Line', 'Tag', tag);
    if isempty(h) || ~isvalid(h)
        h = plot(ax, xPlot, yPlot, 'o', ...
            'LineStyle', 'none', ...
            'MarkerFaceColor', color, ...
            'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'MarkerSize', 5, ...
            'Tag', tag, ...
            'DisplayName', char(names(i)));
    else
        set(h, 'XData', xPlot, 'YData', yPlot);
        if isprop(h, 'DisplayName') && isempty(h.DisplayName)
            h.DisplayName = char(names(i));
        end
    end
    handles(end+1) = h; %#ok<AGROW>
    labels{end+1} = char(names(i)); %#ok<AGROW>
end
end

function prettyName = prettifyTag(tagName)
prettyName = regexprep(tagName, '([a-z0-9])([A-Z])', '$1 $2');
prettyName = strrep(prettyName, '_', ' ');
prettyName = strtrim(prettyName);
end

function [xOut, yOut] = alignSeries(x, y)
% Trim vectors to common length.
if isempty(x) || isempty(y)
    xOut = [];
    yOut = [];
    return;
end
n = min(numel(x), numel(y));
xOut = x(1:n);
yOut = y(1:n);
xOut = xOut(:);
yOut = yOut(:);
end

function xMeas = getMeasuredPositions(state, meas)
% Determine x-axis coordinates for measured data.

xMeas = meas.position;
xMeas = max(0, min(state.Lp, xMeas));
xMeas = xMeas(:)';
end
