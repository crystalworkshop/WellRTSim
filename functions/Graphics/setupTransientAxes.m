function axesStruct = setupTransientAxes(parentTab, pressureUnitLabel)
% setupTransientAxes Create axes for transient profiles.

    gl = uigridlayout(parentTab, [2, 4]);
    gl.ColumnWidth = {'1x','0.8x','1x','0.9x'};
    gl.RowHeight   = {'1x','1x'};
    gl.Padding     = [10 10 10 10];

    bg = [0.15 0.15 0.15]; fg = [0.8 0.8 0.8]; gridc = [0.4 0.4 0.4];
    titleFS = 11; labelFS = 9; axesFS = 9;

    % Row 1 axes
    ax1 = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, 1, 1, 'Pressure', 'Position [m]', ...
        sprintf('P [%s]', pressureUnitLabel));
    ax2 = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, 1, 2, 'Temperature', 'Position [m]', 'T [°C]');
    ax3 = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, 1, 3, 'Velocities', 'Position [m]', 'u [m/s]');

    axGeom = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, [1 2], 4, 'Geometry', 'Diameter [m]', 'Position [m]');
    axGeom.UserData = struct('lastDp', [], 'lastFractions', []);

    % Row 2 axes
    ax4 = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, 2, 1, 'Saturation Index', 'Position [m]', 'SI [-]');
    ax5 = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, 2, 2, ...
        'Mass Concentration', 'Position [m]', 'c [ppm]');
    ax5.YScale = 'log';
    ax5.XMinorGrid = 'off';
    ax5.YMinorGrid = 'off';
    ax5.YLim = [1e-2 1];
    ax6 = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, 2, 3, 'Flow Rates', 'Position [m]', 'Mass flow [kg/s]');

    axesStruct = struct('ax1',ax1,'ax2',ax2,'ax3',ax3,'ax4',ax4,'ax5',ax5,'ax6',ax6,'axGeom',axGeom);
end

function ax = makeAxis(gl, bg, fg, gridc, axesFS, titleFS, labelFS, row, col, ttl, xl, yl)
    ax = uiaxes(gl); 
    ax.Layout.Row = row; 
    ax.Layout.Column = col;
    ax.Color = bg; ax.XColor = fg; ax.YColor = fg; ax.GridColor = gridc; ax.FontSize = axesFS;
    grid(ax, 'on');
    title(ax, ttl, 'FontSize', titleFS, 'Color', [0.9 0.9 0.9]);
    xlabel(ax, xl, 'FontSize', labelFS, 'Color', fg);
    ylabel(ax, yl, 'FontSize', labelFS, 'Color', fg);
    hold(ax, 'on');
end
