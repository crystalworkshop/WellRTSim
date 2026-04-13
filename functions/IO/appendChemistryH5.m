function state = appendChemistryH5(state)
% appendChemistryH5 Append latest chemistry time-series sample to HDF5.

    if ~state.h5_initialized
        return;
    end
    file = state.h5_file;
    % Determine last filled index
    iFilled = state.stepIndex;
    j = state.h5_chem_idx + 1;
    if j > iFilled
        return; % nothing new to write
    end

    h5write(file, '/chemistry/time_days',  state.tsav(j),         j, 1);
    h5write(file, '/chemistry/iterations', state.ChemIterations(j), j, 1);
    h5write(file, '/chemistry/tol1',       state.ChemTol1(j),      j, 1);
    state.h5_chem_idx = j;
end
