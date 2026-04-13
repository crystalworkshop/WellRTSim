function rv_launchResultView(targetPath, toolDir)
% rv_launchResultView Resolve the target file and open the standalone viewer.

    if nargin < 1
        targetPath = "";
    end
    if nargin < 2 || strlength(string(toolDir)) == 0
        toolDir = fileparts(mfilename('fullpath'));
    end

    projectRoot = fileparts(toolDir);
    resultsFile = resolveResultsFile(targetPath, projectRoot);
    if strlength(resultsFile) == 0
        fprintf('result_view cancelled by user.\n');
        return;
    end

    viewerData = rv_loadResultsH5(resultsFile);
    rv_createResultsViewer(viewerData);
end

function resultsFile = resolveResultsFile(targetPath, projectRoot)
    targetPath = string(targetPath);
    if strlength(strtrim(targetPath)) > 0
        candidate = resolveCandidate(targetPath);
        if strlength(candidate) > 0
            resultsFile = candidate;
            return;
        end
        warning('result_view:MissingFile', ...
            'Specified path did not resolve to a results HDF5 file: %s', targetPath);
    end

    defaultDir = fullfile(projectRoot, 'Simulation');
    if ~isfolder(defaultDir)
        altDir = fullfile(projectRoot, 'Simulations');
        if isfolder(altDir)
            defaultDir = altDir;
        else
            defaultDir = projectRoot;
        end
    end
    [fileName, pathName] = uigetfile({'*.h5', 'HDF5 files (*.h5)'}, ...
        'Select results HDF5', defaultDir);
    if isequal(fileName, 0)
        resultsFile = "";
    else
        resultsFile = string(fullfile(pathName, fileName));
    end
end

function candidate = resolveCandidate(targetPath)
    candidate = "";
    targetPath = string(targetPath);
    if isfolder(targetPath)
        candidate = findLatestResultsFile(char(targetPath));
        return;
    end
    if isfile(targetPath)
        candidate = targetPath;
    end
end

function resultsFile = findLatestResultsFile(simDir)
    resultsFile = "";
    exactFile = fullfile(simDir, 'results.h5');
    listing = dir(fullfile(simDir, 'results*.h5'));
    if isempty(listing)
        if isfile(exactFile)
            resultsFile = string(exactFile);
        end
        return;
    end

    indexedFiles = strings(0, 1);
    indexedIds = zeros(0, 1);
    for k = 1:numel(listing)
        name = string(listing(k).name);
        token = regexp(name, '^results_(\d+)\.h5$', 'tokens', 'once');
        if isempty(token)
            continue;
        end
        indexedFiles(end+1, 1) = string(fullfile(listing(k).folder, listing(k).name)); %#ok<AGROW>
        indexedIds(end+1, 1) = str2double(token{1}); %#ok<AGROW>
    end

    if ~isempty(indexedFiles)
        [~, idx] = max(indexedIds);
        resultsFile = indexedFiles(idx);
        return;
    end

    if isfile(exactFile)
        resultsFile = string(exactFile);
        return;
    end

    [~, idx] = max([listing.datenum]);
    resultsFile = string(fullfile(listing(idx).folder, listing(idx).name));
end
