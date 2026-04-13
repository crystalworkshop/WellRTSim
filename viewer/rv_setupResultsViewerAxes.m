function axesStruct = rv_setupResultsViewerAxes(parentContainer)
% rv_setupResultsViewerAxes Create the UI axes layout for the results viewer.

    bg = [0.15 0.15 0.15];
    fg = [0.85 0.85 0.85];
    gridc = [0.4 0.4 0.4];
    titlec = [0.95 0.95 0.95];
    titleFS = 11;
    labelFS = 9;
    axesFS = 9;

    outer = uigridlayout(parentContainer, [1, 3]);
    outer.ColumnWidth = {'4x', '1.1x', '1.4x'};
    outer.RowHeight = {'1x'};
    outer.Padding = [0 0 0 0];
    outer.BackgroundColor = bg;

    left = uigridlayout(outer, [2, 3]);
    left.Layout.Row = 1;
    left.Layout.Column = 1;
    left.RowHeight = {'1x', '1x'};
    left.ColumnWidth = {'1x', '0.9x', '1x'};
    left.Padding = [10 10 10 10];
    left.BackgroundColor = bg;

    axesStruct = struct();
    axesStruct.ax1 = makeAxis(left, bg, fg, gridc, axesFS, titleFS, labelFS, 1, 1, ...
        'Pressure', 'Position [m]', 'P', titlec);
    axesStruct.ax2 = makeAxis(left, bg, fg, gridc, axesFS, titleFS, labelFS, 1, 2, ...
        'Temperature', 'Position [m]', 'T [degC]', titlec);
    axesStruct.ax3 = makeAxis(left, bg, fg, gridc, axesFS, titleFS, labelFS, 1, 3, ...
        'Velocities', 'Position [m]', 'u [m/s]', titlec);
    axesStruct.ax4 = makeAxis(left, bg, fg, gridc, axesFS, titleFS, labelFS, 2, 1, ...
        'Saturation Index', 'Position [m]', 'SI [-]', titlec);
    axesStruct.ax5 = makeAxis(left, bg, fg, gridc, axesFS, titleFS, labelFS, 2, 2, ...
        'Mass Concentration', 'Position [m]', 'c [ppm]', titlec);
    axesStruct.ax6 = makeAxis(left, bg, fg, gridc, axesFS, titleFS, labelFS, 2, 3, ...
        'Flow Rates', 'Position [m]', 'Mass flow [kg/s]', titlec);

    geomAx = uiaxes(outer);
    geomAx.Layout.Row = 1;
    geomAx.Layout.Column = 2;
    geomAx.Color = bg;
    geomAx.XColor = fg;
    geomAx.YColor = fg;
    geomAx.GridColor = gridc;
    geomAx.FontSize = axesFS;
    geomAx.BackgroundColor = bg;
    grid(geomAx, 'on');
    title(geomAx, 'Geometry', 'FontSize', titleFS, 'Color', titlec);
    xlabel(geomAx, 'Diameter [m]', 'FontSize', labelFS, 'Color', fg);
    ylabel(geomAx, 'Position [m]', 'FontSize', labelFS, 'Color', fg);
    hold(geomAx, 'on');
    axesStruct.axGeom = geomAx;

    rightPanel = uigridlayout(outer, [4, 2]);
    rightPanel.Layout.Row = 1;
    rightPanel.Layout.Column = 3;
    rightPanel.ColumnWidth = {60, '1x'};
    rightPanel.RowHeight = {'0.9x', '0.6x', '0.9x', 'fit'};
    rightPanel.Padding = [0 0 0 0];
    rightPanel.ColumnSpacing = 0;
    rightPanel.RowSpacing = 0;
    rightPanel.BackgroundColor = bg;

    depthSlider = uislider(rightPanel, ...
        'Orientation', 'vertical', ...
        'Limits', [1 2], ...
        'Value', 1, ...
        'MajorTicksMode', 'manual', ...
        'MinorTicks', []);
    depthSlider.Layout.Row = [1 3];
    depthSlider.Layout.Column = 1;
    depthSlider.Enable = 'off';

    axScaleThickness = uiaxes(rightPanel);
    axScaleThickness.Layout.Row = 1;
    axScaleThickness.Layout.Column = 2;
    axScaleThickness.Color = bg;
    axScaleThickness.XColor = fg;
    axScaleThickness.YColor = fg;
    axScaleThickness.GridColor = gridc;
    axScaleThickness.FontSize = axesFS;
    axScaleThickness.BackgroundColor = bg;
    grid(axScaleThickness, 'on');
    title(axScaleThickness, 'Scale Thickness History', 'FontSize', titleFS, 'Color', titlec);
    xlabel(axScaleThickness, 'Time', 'FontSize', labelFS, 'Color', fg);
    ylabel(axScaleThickness, 'Thickness [m]', 'FontSize', labelFS, 'Color', fg);
    hold(axScaleThickness, 'on');

    axFlow = uiaxes(rightPanel);
    axFlow.Layout.Row = 2;
    axFlow.Layout.Column = 2;
    axFlow.Color = bg;
    axFlow.XColor = fg;
    axFlow.YColor = fg;
    axFlow.GridColor = gridc;
    axFlow.FontSize = axesFS;
    axFlow.BackgroundColor = bg;
    grid(axFlow, 'on');
    title(axFlow, 'Scale Accumulation Rate', 'FontSize', titleFS, 'Color', titlec);
    xlabel(axFlow, 'Time', 'FontSize', labelFS, 'Color', fg);
    ylabel(axFlow, 'Precipitation [kg/s]', 'FontSize', labelFS, 'Color', fg);
    hold(axFlow, 'on');

    axScale = uiaxes(rightPanel);
    axScale.Layout.Row = 3;
    axScale.Layout.Column = 2;
    axScale.Color = bg;
    axScale.XColor = fg;
    axScale.YColor = fg;
    axScale.GridColor = gridc;
    axScale.FontSize = axesFS;
    axScale.BackgroundColor = bg;
    grid(axScale, 'on');
    title(axScale, 'Scale Composition', 'FontSize', titleFS, 'Color', titlec);
    xlabel(axScale, 'Normalized thickness [-]', 'FontSize', labelFS, 'Color', fg);
    ylabel(axScale, 'Mineral fraction [-]', 'FontSize', labelFS, 'Color', fg);
    hold(axScale, 'on');
    axScale.XLim = [0 1];
    axScale.YLim = [0 1];
    axScale.XDir = 'reverse';

    exportButton = uibutton(rightPanel, ...
        'Text', 'Export scale fractions', ...
        'Enable', 'off');
    exportButton.Layout.Row = 4;
    exportButton.Layout.Column = 2;
    exportButton.Tooltip = 'Save scale fractions at selected depth to CSV';

    axesStruct.axScale = axScale;
    axesStruct.axScaleThickness = axScaleThickness;
    axesStruct.axFlow = axFlow;
    axesStruct.depthSlider = depthSlider;
    axesStruct.exportScaleButton = exportButton;
end

function ax = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, row, col, ttl, xl, yl, titlec)
    ax = uiaxes(gl);
    ax.Layout.Row = row;
    ax.Layout.Column = col;
    ax.Color = bg;
    ax.XColor = fg;
    ax.YColor = fg;
    ax.GridColor = gridc;
    ax.FontSize = axesFS;
    ax.BackgroundColor = bg;
    grid(ax, 'on');
    title(ax, ttl, 'FontSize', titleFS, 'Color', titlec);
    xlabel(ax, xl, 'FontSize', labelFS, 'Color', fg);
    ylabel(ax, yl, 'FontSize', labelFS, 'Color', fg);
    hold(ax, 'on');
end
