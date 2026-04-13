function state = appendWellheadH5(state)
% appendWellheadH5 Append latest wellhead time-series sample to HDF5

    if ~state.h5_initialized
        return;
    end
    file = state.h5_file;
    % Determine last filled index
    iFilled = state.stepIndex;
    j = state.h5_well_idx + 1;
    if j > iFilled
        return; % nothing new to write
    end

    % Write one new sample at index j
    h5write(file, '/wellhead/time_days',      state.tsav(j),               j, 1);
    h5write(file, '/wellhead/flow_rate_kgps', state.FlowRate(j),           j, 1);
    pScale = state.pressureUnitScale;
    h5write(file, '/wellhead/pressure_MPa',   state.WHP(j)/pScale,            j, 1);
    h5write(file, '/wellhead/rho_mix',        state.RhoTop(j),             j, 1);
    h5write(file, '/wellhead/rho_gas',    state.RhoGasTop(j),          j, 1);
    h5write(file, '/wellhead/rho_liq',    state.RhoLiqTop(j),          j, 1);
    h5write(file, '/wellhead/quality_mass',   state.Quality(j),            j, 1);
    h5write(file, '/wellhead/u_mix',          state.UMixTop(j),            j, 1);
    h5write(file, '/wellhead/u_gas',          state.UGasTop(j),            j, 1);
    h5write(file, '/wellhead/u_liq',          state.ULiqTop(j),            j, 1);
    h5write(file, '/wellhead/T_C',        state.TTop(j) - 273.15,      j, 1);
    h5write(file, '/wellhead/iterations', state.Iterations(j), j, 1);
    h5write(file, '/wellhead/tol1', state.Tol1(j), j, 1);

    if isfield(state, 'h5_well_chem_paths') && ~isempty(state.h5_well_chem_paths) ...
            && isfield(state, 'CoutPpm') && size(state.CoutPpm, 2) >= j
        chemVals = state.CoutPpm(:, j);
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
