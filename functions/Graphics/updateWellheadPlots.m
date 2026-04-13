function updateWellheadPlots(state)
% updateWellheadPlots Update Well head time-series lines (in-place)
% Expects state.WellAxes with line handles and state vectors:
%   tsav [days], FlowRate [kg/s], WHP [Pa], Quality [-],
%   UMixTop [m/s], UGasTop [m/s], ULiqTop [m/s]
% WHP.csv series are overlaid (pressure/temp/flows) starting at t_adjust+t_transit.

    tfun = @(v) v(:)';
    nUsed = max(1, state.stepIndex);
    idx = 1:nUsed;
    t = tfun(state.tsav(idx));

    sched = state.whpSchedule;
    tDays = (state.t_adjust + state.t_transit + sched.tSec(:)) / 86400;

    flowTotal = sched.flowTotal(:);
    flowBrine = sched.flowBrine(:);
    flowSteam = sched.flowSteam(:);
    pPa = sched.pPa(:);
    tempC = sched.tempC(:);

    set(state.WellAxes.hFlow, 'XData', t, 'YData', tfun(state.FlowRate(idx)));
    set(state.WellAxes.hFlowObsTot, 'XData', tfun(tDays), 'YData', tfun(flowTotal));
    set(state.WellAxes.hFlowObsBrine, 'XData', tfun(tDays), 'YData', tfun(flowBrine));
    set(state.WellAxes.hFlowObsSteam, 'XData', tfun(tDays), 'YData', tfun(flowSteam));

    pScale = state.pressureUnitScale;
    pLabel = char(state.pressureUnitLabel);
    set(state.WellAxes.hPres, 'XData', t, 'YData', tfun(state.WHP(idx)) / pScale);
    set(state.WellAxes.hPresWHP, 'XData', tfun(tDays), 'YData', tfun(pPa) / pScale);
    ylabel(state.WellAxes.axPres, sprintf('P [%s]', pLabel));

    set(state.WellAxes.hSteamCalc, 'XData', t, 'YData', tfun(state.SteamFrac(idx)));
    steamObs = flowSteam(:) ./ flowTotal(:);
    set(state.WellAxes.hSteamObs, 'XData', tfun(tDays), 'YData', tfun(steamObs));

    set(state.WellAxes.hUMix, 'XData', t, 'YData', tfun(state.UMixTop(idx)));
    set(state.WellAxes.hUGas, 'XData', t, 'YData', tfun(state.UGasTop(idx)));
    set(state.WellAxes.hULiq, 'XData', t, 'YData', tfun(state.ULiqTop(idx)));
    set(state.WellAxes.hTemp, 'XData', t, 'YData', tfun(state.TTop(idx) - 273.15));
    set(state.WellAxes.hTempWHP, 'XData', tfun(tDays), 'YData', tfun(tempC));

    chemVals = state.CoutPpm(:, idx);
    nChem = size(chemVals, 1);
    hChem = state.WellAxes.hChem;
    nPlot = min(nChem, numel(hChem));
    if nPlot > 0
        xCell = repmat({t}, nPlot, 1);
        yCell = num2cell(chemVals(1:nPlot, :), 2);
        set(hChem(1:nPlot), {'XData'}, xCell, {'YData'}, yCell);
    end
    if nPlot < numel(hChem)
        emptyCell = repmat({[]}, numel(hChem) - nPlot, 1);
        set(hChem(nPlot+1:end), {'XData'}, emptyCell, {'YData'}, emptyCell);
    end
    ylim(state.WellAxes.axChem, [1e-2 max(1e-1, max(chemVals(:)))]);

    tFinal = state.tfin / 86400;
    if tFinal <= 0
        tFinal = max(t);
    end
    xl = [0, tFinal];
    axList = [state.WellAxes.axFlow, state.WellAxes.axPres, state.WellAxes.axChem, ...
        state.WellAxes.axSteam, state.WellAxes.axVel, state.WellAxes.axTemp];
    axList = axList(isvalid(axList));
    set(axList, 'XLim', xl);
end
