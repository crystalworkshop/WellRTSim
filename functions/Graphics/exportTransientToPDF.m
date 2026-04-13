function exportTransientToPDF(state, outFile)
% Export the transient UIAxes content into a standard figure and save as PDF.
% This avoids UIFigure export limitations and ensures a reliable PDF.

    % Create an invisible standard figure with 2x4 tiled layout
    % Use inches for predictable sizing in PDF and painters for vector output
    f = figure('Visible','off','Color','w','Units','inches','Position',[1 1 12 7.5], ...
               'Renderer','painters');
    tl = tiledlayout(f,2,4,'Padding','compact','TileSpacing','compact');

    % Helper to copy one UIAxes to a standard axes
    function copyAxes(uiAx, titleText)
        ax = nexttile(tl);
        % White background axis for export
        ax.Color = 'w';
        hold(ax,'on'); grid(ax,'on');
        % Copy line and patch children to preserve filled regions
        chL = findall(uiAx,'Type','line');
        chP = findall(uiAx,'Type','patch');
        if ~isempty(chP)
            copyobj(chP, ax);
        end
        if ~isempty(chL)
            newCh = copyobj(chL, ax); %#ok<NASGU>
            try, set(newCh, 'MarkerFaceColor','none'); end
        end
        % Labels and title
        try, ax.XLabel.String = uiAx.XLabel.String; end
        try, ax.YLabel.String = uiAx.YLabel.String; end
        try, ax.Title.String  = uiAx.Title.String; end
        if nargin>1 && ~isempty(titleText)
            ax.Title.String = titleText;
        end
        % Styling: smaller fonts, black text, white background
        ax.XColor = [0 0 0];
        ax.YColor = [0 0 0];
        ax.GridColor = 0.75*[1 1 1];
        ax.Title.Color = [0 0 0];
        ax.XLabel.Color = [0 0 0];
        ax.YLabel.Color = [0 0 0];
        % Thinner axes lines and shorter ticks
        ax.LineWidth = 0.6;
        try, ax.GridAlpha = 0.2; end
        try, ax.MinorGridAlpha = 0.15; end
        try, ax.TickLength = [0.008 0.008]; end
        % Smaller fonts for PDF only and disable multipliers/bold
        ax.TitleFontSizeMultiplier = 1.0;
        ax.LabelFontSizeMultiplier = 1.0;
        ax.FontSize = 6;                 % ticks
        ax.XLabel.FontSize = 6;
        ax.YLabel.FontSize = 6;
        ax.Title.FontSize  = 7;          % smaller tile headers
        try, ax.Title.FontWeight = 'normal'; end
        try, ax.Title.FontName = 'Helvetica'; end
        try, ax.FontName = 'Helvetica'; end
        hold(ax,'off');
    end

    % Copy each transient axis
    copyAxes(state.TranAxes.ax1);
    copyAxes(state.TranAxes.ax2);
    copyAxes(state.TranAxes.ax3);
    copyAxes(state.TranAxes.ax4);
    copyAxes(state.TranAxes.ax5);
    copyAxes(state.TranAxes.ax6);

    % If geometry axis exists, copy it as the last tile
    copyAxes(state.TranAxes.axGeom);

    % Export to vector PDF with the same font settings and white background
    exportgraphics(f, outFile, 'ContentType','vector', 'BackgroundColor','white');
    close(f);
end
