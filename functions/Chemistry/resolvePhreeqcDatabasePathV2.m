function dbPath = resolvePhreeqcDatabasePathV2(dbPath)
% resolvePhreeqcDatabasePathV2 Resolve a PHREEQC database file on disk.

dbPath = strtrim(string(dbPath));
if strlength(dbPath) == 0
    error('resolvePhreeqcDatabasePathV2:EmptyDatabasePath', ...
        'PHREEQC database path must not be empty.');
end

candidate = char(dbPath);
if should_search_matlab_path(candidate)
    if ~isfile(candidate)
        pathDirs = strsplit(path, pathsep);
        for i = 1:numel(pathDirs)
            probe = fullfile(pathDirs{i}, candidate);
            if isfile(probe)
                candidate = probe;
                break
            end
        end
    end
end

if ~isfile(candidate)
    error('resolvePhreeqcDatabasePathV2:DatabaseNotFound', ...
        ['PHREEQC database "%s" was not found. Run the external ' ...
        'PhreeqcMatlab startup.m first or set an absolute path in ' ...
        'chemistry.md.'], char(dbPath));
end

dbPath = string(char(java.io.File(candidate).getCanonicalPath()));
end

function tf = should_search_matlab_path(candidate)
tf = isempty(regexp(candidate, '^[A-Za-z]:[\\/]', 'once')) && ...
    ~startsWith(candidate, filesep) && ...
    ~contains(candidate, '\') && ...
    ~contains(candidate, '/');
end
