function axesStruct = setupWellheadAxes(parentTab, pressureUnitLabel, chemLabels, gasMask)
% setupWellheadAxes Create time-series axes for wellhead diagnostics
% Axes: Flow rate, Pressure, Outlet chemistry, Steam fraction, Velocities (mix/g/l)

    chemLabels = string(chemLabels(:));
    gasMask = gasMask(:);

    gl = uigridlayout(parentTab, [2, 3]);
    gl.ColumnWidth = {'1x','1x','1x'};
    gl.RowHeight   = {'1x','1x'};
    gl.Padding = [10 10 10 10];

    bg = [0.15 0.15 0.15]; fg = [0.85 0.85 0.85]; gridc = [0.4 0.4 0.4];
    titlec = [0.95 0.95 0.95]; titleFS = 11; labelFS = 9; axesFS = 9;

    % Helper to configure UIAxes
    function ax = makeAx(row,col,ttl,xl,yl)
        ax = uiaxes(gl); ax.Layout.Row = row; ax.Layout.Column = col;
        ax.Color = bg; ax.XColor = fg; ax.YColor = fg; ax.GridColor = gridc; ax.FontSize = axesFS;
        grid(ax,'on'); title(ax, ttl, 'FontSize', titleFS, 'Color', titlec);
        xlabel(ax, xl, 'FontSize', labelFS, 'Color', fg);
        ylabel(ax, yl, 'FontSize', labelFS, 'Color', fg);
        hold(ax,'on');
    end

    axFlow = makeAx(1,1,'Flow Rate','Time [days]','m_dot [kg/s]');
    axPres = makeAx(1,2,'Wellhead Pressure','Time [days]',sprintf('P [%s]', pressureUnitLabel));
    axChem = makeAx(1,3,'Outlet Mass Concentration','Time [days]','c [ppm]');
    axChem.YScale = 'log';
    axChem.XMinorGrid = 'off';
    axChem.YMinorGrid = 'off';
    axChem.YLim = [1e-2 1];
    axSteam = makeAx(2,1,'Steam Fraction','Time [days]','x_{steam} [-]');
    axVel  = makeAx(2,2,'Velocities','Time [days]','u [m/s]');
    % Temperature axis in (2,3)
    axTemp = makeAx(2,3,'Temperature','Time [days]','T [°C]');

    % Create line placeholders
    t0 = 0; y0 = 0;
    markerSize = 3;
    hFlowObsTot = plot(axFlow, t0, y0, 'o', 'LineStyle', 'none', 'MarkerSize', markerSize, ...
        'MarkerFaceColor', [0.9 0.7 0.2], 'MarkerEdgeColor', [0.9 0.7 0.2], ...
        'DisplayName', 'WHP total');
    hFlowObsBrine = plot(axFlow, t0, y0, 's', 'LineStyle', 'none', 'MarkerSize', markerSize, ...
        'MarkerFaceColor', [0.3 0.9 0.9], 'MarkerEdgeColor', [0.3 0.9 0.9], ...
        'DisplayName', 'WHP brine');
    hFlowObsSteam = plot(axFlow, t0, y0, '^', 'LineStyle', 'none', 'MarkerSize', markerSize, ...
        'MarkerFaceColor', [1.0 0.5 0.2], 'MarkerEdgeColor', [1.0 0.5 0.2], ...
        'DisplayName', 'WHP steam');
    hFlow = plot(axFlow, t0, y0, '-', 'LineWidth', 1.2, 'Color', [0.9 0.9 0.2], ...
        'DisplayName', 'Sim total');
    hPresWHP = plot(axPres, t0, y0, 'o', 'LineStyle', 'none', 'MarkerSize', markerSize, ...
        'MarkerFaceColor', [0.9 0.7 0.2], 'MarkerEdgeColor', [0.9 0.7 0.2], ...
        'DisplayName', 'WHP');
    hPres = plot(axPres, t0, y0, '-', 'LineWidth', 1.2, 'Color', [0.2 0.8 0.9], ...
        'DisplayName', 'Simulated');
    hSteamObs = plot(axSteam, t0, y0, 'o', 'LineStyle', 'none', 'MarkerSize', markerSize, ...
        'MarkerFaceColor', [0.9 0.7 0.2], 'MarkerEdgeColor', [0.9 0.7 0.2], ...
        'DisplayName', 'WHP');
    hSteamCalc = plot(axSteam, t0, y0, '-', 'LineWidth', 1.2, 'Color', [0.9 0.6 0.2], ...
        'DisplayName', 'Calc');
    hUMix = plot(axVel,  t0, y0, '-', 'LineWidth', 1.2, 'Color', [0.2 0.9 0.2], ...
        'DisplayName', 'u_mix');
    hUGas = plot(axVel,  t0, y0, '--','LineWidth', 1.0, 'Color', [1.0 0.3 0.3], ...
        'DisplayName', 'u_g');
    hULiq = plot(axVel,  t0, y0, '--','LineWidth', 1.0, 'Color', [0.3 0.9 0.9], ...
        'DisplayName', 'u_l');
    hTempWHP = plot(axTemp, t0, y0, 'o', 'LineStyle', 'none', 'MarkerSize', markerSize, ...
        'MarkerFaceColor', [0.9 0.7 0.2], 'MarkerEdgeColor', [0.9 0.7 0.2], ...
        'DisplayName', 'WHP');
    hTemp = plot(axTemp, t0, y0, '-', 'LineWidth', 1.2, 'Color', [1.0 0.5 0.2], ...
        'DisplayName', 'Simulated');

    legFlow = legend(axFlow, [hFlow,hFlowObsTot,hFlowObsBrine,hFlowObsSteam], ...
        {'Sim total','WHP total','WHP brine','WHP steam'}, 'TextColor', fg, 'Location','best');
    legSteam = legend(axSteam, [hSteamCalc,hSteamObs], {'Calc','WHP'}, 'TextColor', fg, 'Location','best');
    legVel = legend(axVel, [hUMix,hUGas,hULiq], {'u_{mix}','u_g','u_l'}, 'TextColor', fg, 'Location','best');
    if isprop(legFlow, 'ItemTokenSize'), legFlow.ItemTokenSize = [8 8]; end
    if isprop(legSteam, 'ItemTokenSize'), legSteam.ItemTokenSize = [8 8]; end
    if isprop(legVel, 'ItemTokenSize'), legVel.ItemTokenSize = [8 8]; end

    hChem = gobjects(0, 1);
    if ~isempty(chemLabels)
        colors = chemPlotPalette(numel(chemLabels));
        hChem = gobjects(numel(chemLabels), 1);
        yChem0 = NaN;
        for i = 1:numel(chemLabels)
            lineStyle = '-';
            if gasMask(i)
                lineStyle = '--';
            end
            hChem(i) = semilogy(axChem, t0, yChem0, lineStyle, 'LineWidth', 1.2, ...
                'Color', colors(i, :), 'DisplayName', char(chemLabels(i)));
        end
        legChem = legend(axChem, hChem, cellstr(chemLabels), 'TextColor', fg, 'Location', 'best');
        if isprop(legChem, 'ItemTokenSize'), legChem.ItemTokenSize = [8 8]; end
    end

    axesStruct = struct(...
        'axFlow',axFlow,'axPres',axPres,'axChem',axChem,'axSteam',axSteam,'axVel',axVel,'axTemp',axTemp, ...
        'hFlow',hFlow,'hFlowObsTot',hFlowObsTot,'hFlowObsBrine',hFlowObsBrine,'hFlowObsSteam',hFlowObsSteam, ...
        'hPres',hPres,'hPresWHP',hPresWHP, ...
        'hSteamCalc',hSteamCalc,'hSteamObs',hSteamObs, ...
        'hUMix',hUMix,'hUGas',hUGas,'hULiq',hULiq,'hTemp',hTemp,'hTempWHP',hTempWHP, ...
        'hChem',hChem);
end
