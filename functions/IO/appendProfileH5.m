function state = appendProfileH5(state, forceSave)
% appendProfileH5 Append a spatial profile snapshot to HDF5/CSV on schedule.
%
% The save cadence is controlled by state.profile_save_interval (seconds).

    if nargin < 2
        forceSave = false;
    end
    if ~state.h5_initialized
        return;
    end

    interval = max(state.profile_save_interval, 1e-9);
    if ~forceSave && (state.tt + 1e-9) < state.next_profile_save_time
        return;
    end

    file = state.h5_file;
    n = state.n;
    j = state.h5_prof_idx + 1;

    % Time in days
    t_days = state.tt/3600/24;
    % Record time in the profiles timeline
    try, h5write(file, '/profiles/time_days', t_days, j, 1); catch, end

    % Extract current fields
    P = state.Y(1, :);
    H = state.Y(2, :);
    U = state.Y(3, :);

    pScale = state.pressureUnitScale;
    pLabel = char(state.pressureUnitLabel);
    P_out = P / pScale;
    Ycurr = [P; H; U];

    % Cache time in days for file naming consistency with legacy folders
    tDaysStr = sprintf('%.6f', t_days);

    % Compute per-cell properties akin to plotResultsOnAxes
    [alpha_g, alpha_l, rho_g, rho_l, ~, ~, ~, T] = calculatePhaseProperties(P, H, state);
    rho_mix      = alpha_g .* rho_g + alpha_l .* rho_l;
    numQuality   = alpha_g .* rho_g;
    denQuality   = numQuality + (1 - alpha_g) .* rho_l;
    quality_mass = zeros(size(alpha_g));
    validQuality = denQuality > 0;
    quality_mass(validQuality) = numQuality(validQuality) ./ denQuality(validQuality);
    thetaNodes = state.gravityThetaNode(:).';
    [u_g, u_l] = calculatePhaseVelocities(U, alpha_g, rho_g, rho_l, state.Dp(:).', T, thetaNodes, state);
    state = refreshHydroFluxCache(state, Ycurr);
    qGas = state.Q_v(1:n);
    qLiq = state.Q_l(1:n);
    qMix = qGas + qLiq;
    qGasFace = state.Q_v_face(1:n+1);
    qLiqFace = state.Q_l_face(1:n+1);

    % Chemistry overlays
    nMinerals = 0;
    mineralNames = strings(0);
    siMatrix = [];
    scaleFrac = [];
    currThickness = [];
    elementProfiles = [];
    elementSymbols = strings(0);

    mineralNames = string(state.chem.mineralNames(:));
    nMinerals = numel(mineralNames);
    currThickness = zeros(nMinerals, n);
    if ~isempty(state.chem.mineralThickness)
        rr = min(nMinerals, size(state.chem.mineralThickness, 1));
        cc = min(n, size(state.chem.mineralThickness, 2));
        currThickness(1:rr, 1:cc) = state.chem.mineralThickness(1:rr, 1:cc);
    end

    siNames = string(state.chem.mineralNames(:));
    siSource = state.chem.saturationIndices;
    siMatrix = zeros(nMinerals, n);
    [hasSiRow, siRowIdx] = ismember(mineralNames, siNames);
    validSi = hasSiRow & siRowIdx <= size(siSource, 1);
    if any(validSi)
        cc = min(n, size(siSource, 2));
        siMatrix(validSi, 1:cc) = siSource(siRowIdx(validSi), 1:cc);
    end

    prevMatrix = state.profileHistory.mineralThickness;
    if isempty(prevMatrix) || size(prevMatrix, 1) ~= nMinerals || size(prevMatrix, 2) ~= n
        prevAligned = zeros(nMinerals, n);
        if ~isempty(prevMatrix)
            rr = min(nMinerals, size(prevMatrix, 1));
            cc = min(n, size(prevMatrix, 2));
            prevAligned(1:rr, 1:cc) = prevMatrix(1:rr, 1:cc);
        end
        prevMatrix = prevAligned;
    end
    deltaThickness = max(currThickness - prevMatrix, 0);
    totalDelta = sum(deltaThickness, 1);
    scaleFrac = zeros(size(deltaThickness));
    hasDelta = totalDelta > 1e-12;
    if any(hasDelta)
        scaleFrac(:, hasDelta) = deltaThickness(:, hasDelta) ./ totalDelta(hasDelta);
    end
    if any(~hasDelta)
        totalCurrent = sum(currThickness(:, ~hasDelta), 1);
        fallbackMask = totalCurrent > 1e-12;
        if any(fallbackMask)
            idx = find(~hasDelta);
            idx = idx(fallbackMask);
            scaleFrac(:, idx) = currThickness(:, idx) ./ totalCurrent(fallbackMask);
        end
    end
    state.profileHistory.mineralThickness = currThickness;

    elementSymbols = string(state.chem.elementSymbols(:));
    elemIdx = state.chem.elementSpeciesIndex(:);
    nElements = numel(elementSymbols);
    if nElements > 0
        elementProfiles = state.C(elemIdx, 1:n);
    else
        elementProfiles = zeros(0, n);
    end

    % Create an HDF5 subgroup under /profiles with this time as the name
    grpName = sprintf('/profiles/%.6f', t_days);
    % Fixed-size datasets under this group
    h5create(file, [grpName '/P'],       [n], 'Datatype','double');   h5write(file, [grpName '/P'],       P_out.');
    try
        h5writeatt(file, [grpName '/P'], 'units', pLabel);
    catch
    end
    h5create(file, [grpName '/H'],       [n], 'Datatype','double');   h5write(file, [grpName '/H'],       H.');
    try
        h5writeatt(file, [grpName '/H'], 'units', 'J/kg');
        h5writeatt(file, [grpName '/H'], 'description', 'mixture enthalpy');
    catch
    end
    h5create(file, [grpName '/T_C'],     [n], 'Datatype','double');   h5write(file, [grpName '/T_C'],     (T - 273.15).');
    h5create(file, [grpName '/u_mix'],   [n], 'Datatype','double');   h5write(file, [grpName '/u_mix'],   U.');
    h5create(file, [grpName '/u_gas'],   [n], 'Datatype','double');   h5write(file, [grpName '/u_gas'],   u_g.');
    h5create(file, [grpName '/u_liq'],   [n], 'Datatype','double');   h5write(file, [grpName '/u_liq'],   u_l.');
    h5create(file, [grpName '/rho_mix'], [n], 'Datatype','double');   h5write(file, [grpName '/rho_mix'], rho_mix.');
    h5create(file, [grpName '/rho_gas'], [n], 'Datatype','double');   h5write(file, [grpName '/rho_gas'], rho_g.');
    h5create(file, [grpName '/rho_liq'], [n], 'Datatype','double');   h5write(file, [grpName '/rho_liq'], rho_l.');
    h5create(file, [grpName '/Q_mass'],  [n], 'Datatype','double');   h5write(file, [grpName '/Q_mass'],  state.Q_mass.');
    h5create(file, [grpName '/Q_v'],     [n], 'Datatype','double');   h5write(file, [grpName '/Q_v'],     qGas(:).');
    h5create(file, [grpName '/Q_l'],     [n], 'Datatype','double');   h5write(file, [grpName '/Q_l'],     qLiq(:).');
    h5create(file, [grpName '/Q_mix'],   [n], 'Datatype','double');   h5write(file, [grpName '/Q_mix'],   qMix(:).');
    h5create(file, [grpName '/Q_v_face'], [n+1], 'Datatype','double'); h5write(file, [grpName '/Q_v_face'], qGasFace(:).');
    h5create(file, [grpName '/Q_l_face'], [n+1], 'Datatype','double'); h5write(file, [grpName '/Q_l_face'], qLiqFace(:).');
    h5create(file, [grpName '/alpha_g'], [n], 'Datatype','double');   h5write(file, [grpName '/alpha_g'], alpha_g.');
    h5create(file, [grpName '/quality_mass'], [n], 'Datatype','double');
    h5write(file, [grpName '/quality_mass'], quality_mass.');

    % Mineral saturation indices and scale composition
    for idx = 1:nMinerals
        mineralName = char(mineralNames(idx));
        datasetBase = matlab.lang.makeValidName(mineralName);
        siPath = sprintf('%s/SI_%s', grpName, datasetBase);
        try
            h5create(file, siPath, [n], 'Datatype','double');
        catch
        end
        try
            h5write(file, siPath, siMatrix(idx, 1:n).');
        catch
        end
        fracPath = sprintf('%s/scaleFrac_%s', grpName, datasetBase);
        try
            h5create(file, fracPath, [n], 'Datatype','double');
        catch
        end
        try
            h5write(file, fracPath, scaleFrac(idx, 1:n).');
        catch
        end
    end

    % Element profiles (mass fraction / tracked species)
    for idx = 1:numel(elementSymbols)
        elemName = char(elementSymbols(idx));
        datasetBase = matlab.lang.makeValidName(['element_' elemName]);
        elemPath = sprintf('%s/%s', grpName, datasetBase);
        try
            h5create(file, elemPath, [n], 'Datatype','double');
        catch
        end
        try
            h5write(file, elemPath, elementProfiles(idx, 1:n).');
        catch
        end
    end

    % Gas partition coefficients (K = c_g/c_l)
    K = state.chem.partitionCoefficients;
    entries = state.chem.species;
    types = arrayfun(@(e) string(e.type), entries);
    isGas = strcmpi(types, "g") | strcmpi(types, "gas");
    gasEntries = entries(isGas);
    chemNames = string(state.chemNames(:));
    gasNames = string({gasEntries.name});
    [hasCompIdx, compIdx] = ismember(upper(gasNames), upper(chemNames));
    gasNames = gasNames(hasCompIdx);
    compIdx = compIdx(hasCompIdx);
    for idx = 1:numel(gasNames)
        gasName = gasNames(idx);
        datasetBase = matlab.lang.makeValidName(['K_' char(gasName)]);
        kPath = sprintf('%s/%s', grpName, datasetBase);
        try
            h5create(file, kPath, [n], 'Datatype','double');
        catch
        end
        try
            kVals = K(compIdx, 1:n);
            h5write(file, kPath, kVals(:).');
        catch
        end
        try
            h5writeatt(file, kPath, 'units', '-');
            h5writeatt(file, kPath, 'species_name', char(gasName));
            h5writeatt(file, kPath, 'description', 'gas partition coefficient c_g/c_l');
        catch
        end
    end

    % Geometry
    try
        h5create(file, [grpName '/Dp'], [n], 'Datatype','double'); h5write(file, [grpName '/Dp'], state.Dp.');
    catch
    end

    state.h5_prof_idx = j;

    % Also create a per-time folder with a CSV snapshot
    if state.save_csv
        try
            profRoot = fullfile(state.SimDir, 'profiles');
            if ~exist(profRoot,'dir'), mkdir(profRoot); end

            pVar = matlab.lang.makeValidName(['P_' pLabel]);
            Ttbl = table( ...
                state.x(:), P_out(:), (T(:)-273.15), U(:), u_g(:), u_l(:), ...
                rho_mix(:), rho_g(:), rho_l(:), ...
                'VariableNames', {'x_m', pVar, 'T_C', 'u_mix', 'u_gas', 'u_liq', 'rho_mix', 'rho_gas', 'rho_liq'});
            Ttbl.Q_mass = state.Q_mass(:);
            Ttbl.Q_v = qGas(:);
            Ttbl.Q_l = qLiq(:);
            Ttbl.Q_mix = qMix(:);
            Ttbl.alpha_g = alpha_g(:);
            Ttbl.quality_mass = quality_mass(:);
            Ttbl.Dp = state.Dp(:);
            Ttbl.H_Jkg = H(:);
            for idx = 1:nMinerals
                mineralName = char(mineralNames(idx));
                siVar = matlab.lang.makeValidName(['SI_' mineralName]);
                fracVar = matlab.lang.makeValidName(['scaleFrac_' mineralName]);
                Ttbl.(siVar) = siMatrix(idx, 1:n).';
                Ttbl.(fracVar) = scaleFrac(idx, 1:n).';
            end
            for idx = 1:numel(elementSymbols)
                elemName = char(elementSymbols(idx));
                elemVar = matlab.lang.makeValidName(['element_' elemName]);
                Ttbl.(elemVar) = elementProfiles(idx, 1:n).';
            end
            K = state.chem.partitionCoefficients;
            for idx = 1:numel(gasNames)
                gasName = char(gasNames(idx));
                kVar = matlab.lang.makeValidName(['K_' gasName]);
                kVals = K(compIdx(idx), 1:n);
                Ttbl.(kVar) = kVals(:);
            end

            outCsv = fullfile(profRoot, ['profile_' tDaysStr 'd.csv']);
            writetable(Ttbl, outCsv);
        catch ME
            warning('Profile folder write failed: %s', ME.message);
        end
    end

    state.h5_prof_idx = j;
    state.last_profile_save_time = state.tt;
    nextTime = state.next_profile_save_time;
    if forceSave && (state.tt + 1e-9) < nextTime
        nextTime = state.tt;
    end
    while (state.tt + 1e-9) >= nextTime
        nextTime = nextTime + interval;
    end
    state.next_profile_save_time = nextTime;
end
