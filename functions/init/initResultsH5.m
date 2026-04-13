function state = initResultsH5(state)
% initResultsH5 Initialize/prepare an HDF5 results file in state.SimDir
% - Creates extendable datasets for wellhead time series and profile snapshots

    if ~exist(state.SimDir,'dir')
        try, mkdir(state.SimDir); end
    end

    outFile = getNextIndexedFile(state.SimDir, resolveResultsBaseName(state), 'h5', 3);
    state.h5_file = outFile;
    state.h5_initialized = false;
    state.h5_well_idx = 0;
    state.h5_chem_idx = 0;
    state.h5_prof_idx = 0;

    % If file exists (should be unique), delete to start clean
    if exist(outFile,'file') == 2
        try, delete(outFile); catch, end
    end

    n = numel(state.x);

    % Meta datasets (fixed size)
    h5create(outFile, '/meta/x',    [n]);
    h5create(outFile, '/meta/Dp0',  [n]);
    h5write(outFile,  '/meta/x',   state.x(:));
    h5write(outFile, '/meta/Dp0', state.Dp0(:));
    % Meta units (no per-variable descriptions)
    try
        h5writeatt(outFile, '/meta/x',   'units', 'm');
        h5writeatt(outFile, '/meta/Dp0', 'units', 'm');
    catch
    end

    % Wellhead time series (extendable second dimension)
    chunkT = 1024;
    h5create(outFile, '/wellhead/time_days',      [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/flow_rate_kgps', [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/pressure_MPa',   [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/rho_mix',        [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/rho_gas',        [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/rho_liq',        [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/quality_mass',   [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/u_mix',          [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/u_gas',          [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/u_liq',          [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/T_C',            [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/iterations',     [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/wellhead/tol1',           [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    state.h5_well_chem_paths = strings(0, 1);
    % Wellhead units (no per-variable descriptions)
    try
        h5writeatt(outFile, '/wellhead/time_days',      'units', 'days');
        h5writeatt(outFile, '/wellhead/flow_rate_kgps', 'units', 'kg/s');
        h5writeatt(outFile, '/wellhead/pressure_MPa',   'units', state.pressureUnitLabel);
        h5writeatt(outFile, '/wellhead/rho_mix',        'units', 'kg/m^3');
        h5writeatt(outFile, '/wellhead/rho_gas',        'units', 'kg/m^3');
        h5writeatt(outFile, '/wellhead/rho_liq',        'units', 'kg/m^3');
        h5writeatt(outFile, '/wellhead/quality_mass',   'units', '-');
        h5writeatt(outFile, '/wellhead/u_mix',          'units', 'm/s');
        h5writeatt(outFile, '/wellhead/u_gas',          'units', 'm/s');
        h5writeatt(outFile, '/wellhead/u_liq',          'units', 'm/s');
        h5writeatt(outFile, '/wellhead/T_C',            'units', 'degC');
        h5writeatt(outFile, '/wellhead/iterations',     'units', 'count');
        h5writeatt(outFile, '/wellhead/tol1',           'units', '-');
    catch
    end
    try
        state.h5_well_chem_paths = initWellheadChemistryDatasets(outFile, state, chunkT);
    catch
        state.h5_well_chem_paths = strings(0, 1);
    end

    % Chemistry time series (extendable)
    h5create(outFile, '/chemistry/time_days',  [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/chemistry/iterations', [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    h5create(outFile, '/chemistry/tol1',       [Inf], 'ChunkSize', [chunkT], 'Datatype','double');
    try
        h5writeatt(outFile, '/chemistry/time_days',  'units', 'days');
        h5writeatt(outFile, '/chemistry/iterations', 'units', 'count');
        h5writeatt(outFile, '/chemistry/tol1',       'units', '-');
    catch
    end

    % Profiles: only a time list; each snapshot goes in its own subgroup
    h5create(outFile, '/profiles/time_days', [Inf], 'ChunkSize', [256], 'Datatype','double');
    try, h5writeatt(outFile, '/profiles/time_days', 'units', 'days'); catch, end

    % Persist preferred time unit metadata for post-processing
    try
        tunit = char(state.tunit);
        h5writeatt(outFile, '/', 'time_unit', tunit);
        h5writeatt(outFile, '/wellhead/time_days', 'tunit', tunit);
        h5writeatt(outFile, '/chemistry/time_days', 'tunit', tunit);
        h5writeatt(outFile, '/profiles/time_days', 'tunit', tunit);
        h5writeatt(outFile, '/', 'time_unit_seconds', state.timeUnitSeconds);
        h5writeatt(outFile, '/', 'pressure_unit', state.pressureUnitLabel);
        h5writeatt(outFile, '/', 'pressure_unit_scale', state.pressureUnitScale);
    catch
    end

    % Store input files for reproducibility
    try
        writeTextFileToH5(outFile, '/inputs/params_md', state.paramsFile);
        writeTextFileToH5(outFile, '/inputs/chemistry_md', getChemistrySetupPath(state));
        chemTemplatePath = resolveChemistryTemplatePath(state);
        if strlength(chemTemplatePath) > 0
            writeTextFileToH5(outFile, '/inputs/chemistry_pht', char(chemTemplatePath));
        end
    catch
    end

    state.h5_initialized = true;
    % Single text field inside HDF5 with variable descriptions
    try
        pLabel = char(state.pressureUnitLabel);
        desc = [ ...
            "3pointWellbore HDF5 variables summary" newline newline ...
            "[meta]" newline ...
            "  /meta/x        : position [m]" newline ...
            "  /meta/Dp0      : initial inner diameter [m]" newline newline ...
            "[wellhead] (time series)" newline ...
            "  /wellhead/time_days      : time [days]" newline ...
            "  /wellhead/flow_rate_kgps : mass flow rate [kg/s]" newline ...
            "  /wellhead/pressure_MPa   : wellhead pressure [" + string(pLabel) + "]" newline ...
            "  /wellhead/rho_mix        : mixture density [kg/m^3]" newline ...
            "  /wellhead/rho_gas        : gas density [kg/m^3]" newline ...
            "  /wellhead/rho_liq        : liquid density [kg/m^3]" newline ...
            "  /wellhead/quality_mass   : mass quality [-]" newline ...
            "  /wellhead/u_mix          : mixture velocity [m/s]" newline ...
            "  /wellhead/u_gas          : gas velocity [m/s]" newline ...
            "  /wellhead/u_liq          : liquid velocity [m/s]" newline ...
            "  /wellhead/T_C            : temperature [degC]" newline ...
            "  /wellhead/iterations     : hydrodynamics iterations [count]" newline ...
            "  /wellhead/tol1           : hydrodynamics solver tol1 [-]" newline ...
            "  /wellhead/chemistry_ppm_comp_<i>_<name> : outlet component concentration [ppm]" newline ...
            "    attributes: label, species_name, phase, source_index" newline newline ...
            "[chemistry] (time series)" newline ...
            "  /chemistry/time_days     : time [days]" newline ...
            "  /chemistry/iterations    : chemistry iterations [count]" newline ...
            "  /chemistry/tol1          : chemistry solver tol1 [-]" newline newline ...
            "[profiles] (each snapshot under /profiles/<time_days>/...)" newline ...
            "  P [" + string(pLabel) + "], H [J/kg], T_C [degC], u_mix [m/s], u_gas [m/s], u_liq [m/s]," newline ...
            "  rho_mix [kg/m^3], rho_gas [kg/m^3], rho_liq [kg/m^3], Q_mass [kg/(m^3 s)]," newline ...
            "  Q_v/Q_l/Q_mix [kg/s], Q_v_face/Q_l_face [kg/s], Dp [m]," newline ...
            "  element_<name> : tracked chemistry profile" newline ...
            "  K_<gas> : gas partition coefficient (c_g/c_l) [-]" newline ...
            "  H_Jkg : mixture enthalpy in profile CSV snapshots [J/kg]" newline newline ...
            "[inputs]" newline ...
            "  /inputs/params_md : params.md text (UTF-8 string)" newline ...
            "  /inputs/chemistry_md : chemistry setup markdown text (UTF-8 string)" newline ...
            "  /inputs/chemistry_pht : chemistry template text (UTF-8 string)" newline];
        bytes = uint8(desc);
        h5create(outFile, '/description', numel(bytes), 'Datatype','uint8');
        h5write(outFile,  '/description', bytes(:));
        try, h5writeatt(outFile, '/description', 'encoding', 'utf-8'); end
    catch
    end
end

function baseName = resolveResultsBaseName(state)
baseName = strtrim(char(string(state.results_prefix)));
baseName = regexprep(baseName, '\.h5$', '', 'ignorecase');
baseName = regexprep(baseName, '_+$', '');
if isempty(baseName)
    baseName = 'results';
end
end

function writeTextFileToH5(h5file, datasetPath, filePath)
    if exist(filePath, 'file') ~= 2
        return;
    end
    try
        txt = fileread(filePath);
    catch
        return;
    end
    if isempty(txt)
        return;
    end
    % Prefer a string dataset for easy viewing in HDF5 browsers.
    try
        h5create(h5file, datasetPath, 1, 'Datatype', 'string');
        h5write(h5file, datasetPath, string(txt));
        try
            h5writeatt(h5file, datasetPath, 'encoding', 'utf-8');
            h5writeatt(h5file, datasetPath, 'source', filePath);
            h5writeatt(h5file, datasetPath, 'format', 'text');
        catch
        end
        return;
    catch
        % Fallback to raw UTF-8 bytes if string datasets are unsupported.
    end

    bytes = unicode2native(txt, 'UTF-8');
    if isempty(bytes)
        return;
    end
    h5create(h5file, datasetPath, numel(bytes), 'Datatype', 'uint8');
    h5write(h5file, datasetPath, uint8(bytes(:)));
    try
        h5writeatt(h5file, datasetPath, 'encoding', 'utf-8');
        h5writeatt(h5file, datasetPath, 'source', filePath);
        h5writeatt(h5file, datasetPath, 'format', 'text_bytes');
    catch
    end
end

function chemDatasetPaths = initWellheadChemistryDatasets(h5file, state, chunkT)
    labels = string(state.chem.plotLabels(:));
    chemDatasetPaths = strings(numel(labels), 1);
    if isempty(labels)
        return;
    end

    speciesNames = strings(0, 1);
    if isfield(state, 'chemNames') && ~isempty(state.chemNames)
        speciesNames = string(state.chemNames(:));
    end
    gasMask = false(numel(labels), 1);
    if isfield(state, 'chem') && isfield(state.chem, 'gasMask') && ~isempty(state.chem.gasMask)
        gasMask = state.chem.gasMask(:);
        if numel(gasMask) < numel(labels)
            gasMask(end+1:numel(labels), 1) = false;
        else
            gasMask = gasMask(1:numel(labels));
        end
    end

    for idx = 1:numel(labels)
        label = labels(idx);
        datasetBase = matlab.lang.makeValidName(char(label));
        datasetPath = sprintf('/wellhead/chemistry_ppm_comp_%02d_%s', idx, datasetBase);
        chemDatasetPaths(idx) = string(datasetPath);
        h5create(h5file, datasetPath, [Inf], 'ChunkSize', [chunkT], 'Datatype', 'double');
        try
            h5writeatt(h5file, datasetPath, 'units', 'ppm');
            h5writeatt(h5file, datasetPath, 'label', char(label));
            if idx <= numel(speciesNames)
                h5writeatt(h5file, datasetPath, 'species_name', char(speciesNames(idx)));
            end
            if gasMask(idx)
                h5writeatt(h5file, datasetPath, 'phase', 'gas');
            else
                h5writeatt(h5file, datasetPath, 'phase', 'aq');
            end
            h5writeatt(h5file, datasetPath, 'source_index', idx);
        catch
        end
    end
end

function templatePath = resolveChemistryTemplatePath(state)
    templatePath = "";
    chemMdPath = getChemistrySetupPath(state);
    if exist(chemMdPath, 'file') ~= 2
        return;
    end

    try
        cfg = parseChemistrySetup(chemMdPath);
        if isfield(cfg, 'template_path') && strlength(string(cfg.template_path)) > 0
            templatePath = string(cfg.template_path);
        end
    catch
        templatePath = "";
    end

    if strlength(templatePath) > 0
        return;
    end

    fallbackPath = fullfile(state.SimDir, 'chemistry.pht');
    if exist(fallbackPath, 'file') == 2
        templatePath = string(fallbackPath);
    end
end
