function updateWellheadPlots(state)
% updateWellheadPlots Update Well head time-series lines (in-place).
% The simulated history is read back from the HDF5 results file (no in-memory
% time-series arrays are kept). WHP.csv observed series are overlaid.

    tfun = @(v) v(:)';

    % Observed WHP schedule overlays (independent of simulated history)
    sched = state.whpSchedule;
    tDays = (state.t_adjust + state.t_transit + sched.tSec(:)) / 86400;
    flowTotal = sched.flowTotal(:);
    flowBrine = sched.flowBrine(:);
    flowSteam = sched.flowSteam(:);
    pPa = sched.pPa(:);
    tempC = sched.tempC(:);

    pScale = state.pressureUnitScale;
    pLabel = char(state.pressureUnitLabel);

    set(state.WellAxes.hFlowObsTot,   'XData', tfun(tDays), 'YData', tfun(flowTotal));
    set(state.WellAxes.hFlowObsBrine, 'XData', tfun(tDays), 'YData', tfun(flowBrine));
    set(state.WellAxes.hFlowObsSteam, 'XData', tfun(tDays), 'YData', tfun(flowSteam));
    set(state.WellAxes.hPresWHP,      'XData', tfun(tDays), 'YData', tfun(pPa) / pScale);
    steamObs = flowSteam(:) ./ flowTotal(:);
    set(state.WellAxes.hSteamObs,     'XData', tfun(tDays), 'YData', tfun(steamObs));
    set(state.WellAxes.hTempWHP,      'XData', tfun(tDays), 'YData', tfun(tempC));
    ylabel(state.WellAxes.axPres, sprintf('P [%s]', pLabel));

    % Simulated history straight from the HDF5 file
    nw = 0;
    if state.h5_initialized
        nw = state.h5_well_idx;
    end
    if nw >= 1
        file = state.h5_file;
        rd = @(name) reshape(h5read(file, name, 1, nw), 1, []);
        t = rd('/wellhead/time_days');

        set(state.WellAxes.hFlow,      'XData', t, 'YData', rd('/wellhead/flow_rate_kgps'));
        set(state.WellAxes.hPres,      'XData', t, 'YData', rd('/wellhead/pressure_MPa'));
        set(state.WellAxes.hSteamCalc, 'XData', t, 'YData', rd('/wellhead/steam_frac'));
        set(state.WellAxes.hUMix,      'XData', t, 'YData', rd('/wellhead/u_mix'));
        set(state.WellAxes.hUGas,      'XData', t, 'YData', rd('/wellhead/u_gas'));
        set(state.WellAxes.hULiq,      'XData', t, 'YData', rd('/wellhead/u_liq'));
        set(state.WellAxes.hTemp,      'XData', t, 'YData', rd('/wellhead/T_C'));

        hChem = state.WellAxes.hChem;
        chemPaths = strings(0, 1);
        if isfield(state, 'h5_well_chem_paths')
            chemPaths = state.h5_well_chem_paths;
        end
        nPlot = min(numel(chemPaths), numel(hChem));
        chemMax = [];
        for iChem = 1:nPlot
            yChem = rd(char(chemPaths(iChem)));
            set(hChem(iChem), 'XData', t, 'YData', yChem);
            chemMax = max([chemMax, yChem]);
        end
        if nPlot < numel(hChem)
            emptyCell = repmat({[]}, numel(hChem) - nPlot, 1);
            set(hChem(nPlot+1:end), {'XData'}, emptyCell, {'YData'}, emptyCell);
        end
        if ~isempty(chemMax)
            ylim(state.WellAxes.axChem, [1e-2 max(1e-1, chemMax)]);
        end
    end

    tFinal = state.tfin / 86400;
    if tFinal <= 0
        tFinal = max([tDays; 1]);
    end
    xl = [0, tFinal];
    axList = [state.WellAxes.axFlow, state.WellAxes.axPres, state.WellAxes.axChem, ...
        state.WellAxes.axSteam, state.WellAxes.axVel, state.WellAxes.axTemp];
    axList = axList(isvalid(axList));
    set(axList, 'XLim', xl);
end
