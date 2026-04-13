function outFile = getLatestIndexedFile(directory, baseName, ext)
% getLatestIndexedFile Return latest indexed filename or empty if none.

outFile = '';
if nargin < 3
    ext = '';
end

directory = char(directory);
baseName = char(baseName);
ext = char(ext);

if isempty(directory)
    directory = '.';
end
if ~isempty(ext) && ext(1) == '.'
    ext = ext(2:end);
end
if ~ispc && startsWith(directory, '~')
    homeDir = getenv('HOME');
    if ~isempty(homeDir)
        if numel(directory) == 1
            directory = homeDir;
        elseif numel(directory) >= 2 && (directory(2) == '/' || directory(2) == '\\')
            directory = fullfile(homeDir, directory(3:end));
        end
    end
end
if exist(directory, 'dir') ~= 7
    return;
end

pattern = sprintf('%s_*', baseName);
files = dir(fullfile(directory, pattern));
if isempty(files)
    return;
end
maxIdx = -inf;
bestSegCount = inf;
maxName = '';
for i = 1:numel(files)
    if files(i).isdir
        continue;
    end
    name = files(i).name;
    [~, stem, extName] = fileparts(name);
    if isempty(ext)
        if ~isempty(extName)
            continue;
        end
    else
        if isempty(extName) || ~strcmpi(extName(2:end), ext)
            continue;
        end
    end
    [idx, segCount] = parseIndexedStem(stem, baseName);
    if isempty(idx)
        continue;
    end
    if ~isfinite(idx)
        continue;
    end
    if idx > maxIdx || (idx == maxIdx && segCount < bestSegCount)
        maxIdx = idx;
        bestSegCount = segCount;
        maxName = name;
    end
end
if ~isempty(maxName)
    outFile = fullfile(directory, maxName);
end
end

function [idx, segCount] = parseIndexedStem(stem, baseName)
idx = [];
segCount = inf;
if ~startsWith(stem, baseName)
    return;
end
prefixLen = numel(baseName);
if numel(stem) <= prefixLen || stem(prefixLen + 1) ~= '_'
    return;
end
suffix = stem(prefixLen + 2:end);
if isempty(suffix)
    return;
end
parts = strsplit(suffix, '_');
for k = 1:numel(parts)
    if isempty(regexp(parts{k}, '^[0-9]+$', 'once'))
        return;
    end
end
segCount = numel(parts);
idx = str2double(parts{end});
if ~isfinite(idx)
    idx = [];
    segCount = inf;
end
end
