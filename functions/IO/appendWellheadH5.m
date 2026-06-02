function state = appendWellheadH5(state)
% appendWellheadH5 Append the latest wellhead sample (state.wellhead) to HDF5.

    if ~state.h5_initialized
        return;
    end
    file = state.h5_file;
    wh = state.wellhead;
    j = state.h5_well_idx + 1;

    % Write the current sample at the next index
    h5write(file, '/wellhead/time_days',      wh.time_days,        j, 1);
    h5write(file, '/wellhead/flow_rate_kgps', wh.flow_rate,        j, 1);
    pScale = state.pressureUnitScale;
    h5write(file, '/wellhead/pressure_MPa',   wh.WHP/pScale,       j, 1);
    h5write(file, '/wellhead/rho_mix',        wh.rho_mix,          j, 1);
    h5write(file, '/wellhead/rho_gas',        wh.rho_gas,          j, 1);
    h5write(file, '/wellhead/rho_liq',        wh.rho_liq,          j, 1);
    h5write(file, '/wellhead/quality_mass',   wh.quality,          j, 1);
    h5write(file, '/wellhead/steam_frac',     wh.steam_frac,       j, 1);
    h5write(file, '/wellhead/u_mix',          wh.u_mix,            j, 1);
    h5write(file, '/wellhead/u_gas',          wh.u_gas,            j, 1);
    h5write(file, '/wellhead/u_liq',          wh.u_liq,            j, 1);
    h5write(file, '/wellhead/T_C',            wh.T - 273.15,       j, 1);
    h5write(file, '/wellhead/iterations',     wh.iterations,       j, 1);
    h5write(file, '/wellhead/tol1',           wh.tol1,             j, 1);

    if isfield(state, 'h5_well_chem_paths') && ~isempty(state.h5_well_chem_paths) ...
            && ~isempty(wh.CoutPpm)
        chemVals = wh.CoutPpm;
        nChem = min(numel(state.h5_well_chem_paths), numel(chemVals));
        for iChem = 1:nChem
            datasetPath = char(state.h5_well_chem_paths(iChem));
            try
                h5write(file, datasetPath, chemVals(iChem), j, 1);
            catch
            end
        end
    end

    state.h5_well_idx = j;
end
