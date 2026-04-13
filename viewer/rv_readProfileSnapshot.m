function profile = rv_readProfileSnapshot(viewerData, index)
% rv_readProfileSnapshot Load a single profile snapshot from results HDF5.

    arguments
        viewerData (1, 1) struct
        index (1, 1) double {mustBePositive, mustBeInteger}
    end

    if index < 1 || index > numel(viewerData.times)
        error('rv_readProfileSnapshot:IndexOutOfRange', ...
            'Snapshot index %d is out of bounds (1..%d).', ...
            index, numel(viewerData.times));
    end

    file = viewerData.file;
    grp = viewerData.groupNames{index};

    profile = struct();
    profile.timeDays = viewerData.times(index);

    profile.P = readDataset(file, grp, 'P');
    profile.T_C = readDataset(file, grp, 'T_C');
    profile.u_mix = readDataset(file, grp, 'u_mix');
    profile.u_gas = readDataset(file, grp, 'u_gas');
    profile.u_liq = readDataset(file, grp, 'u_liq');
    profile.rho_mix = readDataset(file, grp, 'rho_mix');
    profile.rho_gas = readDataset(file, grp, 'rho_gas');
    profile.rho_liq = readDataset(file, grp, 'rho_liq');

    profile.alpha_g = tryReadDataset(file, grp, 'alpha_g');
    profile.Dp = tryReadDataset(file, grp, 'Dp');
    profile.Q_mass = tryReadDataset(file, grp, 'Q_mass');
    profile.Q_v = tryReadDataset(file, grp, 'Q_v');
    profile.Q_l = tryReadDataset(file, grp, 'Q_l');
    profile.Q_mix = tryReadDataset(file, grp, 'Q_mix');
    profile.Q_v_face = tryReadDataset(file, grp, 'Q_v_face');
    profile.Q_l_face = tryReadDataset(file, grp, 'Q_l_face');
    profile.Qm_l = tryReadDataset(file, grp, 'Qm_l');
    profile.Qm_v = tryReadDataset(file, grp, 'Qm_v');
    profile.Qm_mix = tryReadDataset(file, grp, 'Qm_mix');
    profile.quality_mass = tryReadDataset(file, grp, 'quality_mass');

    profile.SI = cell(0, 1);
    profile.scaleFrac = cell(0, 1);
    profile.elements = cell(0, 1);

    for k = 1:numel(viewerData.siDatasets)
        profile.SI{k} = tryReadDataset(file, grp, viewerData.siDatasets{k}); %#ok<AGROW>
    end
    for k = 1:numel(viewerData.scaleFractionDatasets)
        profile.scaleFrac{k} = tryReadDataset(file, grp, viewerData.scaleFractionDatasets{k}); %#ok<AGROW>
    end
    for k = 1:numel(viewerData.elementDatasets)
        profile.elements{k} = tryReadDataset(file, grp, viewerData.elementDatasets{k}); %#ok<AGROW>
    end
end

function data = readDataset(file, groupName, dataset)
    path = sprintf('%s/%s', groupName, dataset);
    try
        data = h5read(file, path);
    catch ME
        error('rv_readProfileSnapshot:MissingDataset', ...
            'Dataset %s is missing or unreadable (%s).', path, ME.message);
    end
    data = data(:);
end

function data = tryReadDataset(file, groupName, dataset)
    path = sprintf('%s/%s', groupName, dataset);
    try
        data = h5read(file, path);
        data = data(:);
    catch
        data = [];
    end
end
