function state = appendChemistryH5(state)
% appendChemistryH5 Append latest chemistry time-series sample to HDF5.

    if ~state.h5_initialized
        return;
    end
    file = state.h5_file;
    j = state.h5_chem_idx + 1;

    h5write(file, '/chemistry/time_days',  state.wellhead.time_days,    j, 1);
    h5write(file, '/chemistry/iterations', state.chemSample.iterations, j, 1);
    h5write(file, '/chemistry/tol1',       state.chemSample.tol1,       j, 1);
    state.h5_chem_idx = j;
end
