function refreshTransientPlots(stateStruct, Y, titleText)
% Clear transient axes and redraw the current solution snapshot.

axStruct = stateStruct.TranAxes;
fields = fieldnames(axStruct);

for k = 1:numel(fields)
    ax = axStruct.(fields{k});
    cla(ax);
    hold(ax, 'on');
end

plotResultsOnAxes(axStruct, stateStruct, Y, titleText);
end
