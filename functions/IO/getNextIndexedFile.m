function outFile = getNextIndexedFile(directory, baseName, ext, width)
% getNextIndexedFile Return next available indexed filename in directory.
%
% Example: getNextIndexedFile(dir, 'results', 'h5', 3) -> results_001.h5

if nargin < 4 || isempty(width)
    width = 3;
end
if nargin < 3 || isempty(ext)
    ext = '';
end
if isempty(directory)
    directory = '.';
end

baseName = char(baseName);
ext = char(ext);

maxIdx = findMaxIndex(directory, baseName, ext);
nextIdx = maxIdx + 1;
width = max(width, numel(num2str(nextIdx)));

if ~isempty(ext)
    outFile = fullfile(directory, sprintf('%s_%0*d.%s', baseName, width, nextIdx, ext));
else
    outFile = fullfile(directory, sprintf('%s_%0*d', baseName, width, nextIdx));
end
end

function maxIdx = findMaxIndex(directory, baseName, ext)
maxIdx = 0;
maxIdxClean = 0;
maxIdxAny = 0;
if isempty(directory) || exist(directory, 'dir') ~= 7
    return;
end
if isempty(ext)
    pattern = sprintf('%s_*', baseName);
else
    pattern = sprintf('%s_*.%s', baseName, ext);
end
files = dir(fullfile(directory, pattern));
if isempty(files)
    return;
end
for i = 1:numel(files)
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
    if idx > maxIdxAny
        maxIdxAny = idx;
    end
    if segCount == 1 && idx > maxIdxClean
        maxIdxClean = idx;
    end
end
if maxIdxClean > 0
    maxIdx = maxIdxClean;
else
    maxIdx = maxIdxAny;
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
