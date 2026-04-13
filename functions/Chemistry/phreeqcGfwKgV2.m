function mwKg = phreeqcGfwKgV2(formulas, dbPath)
% phreeqcGfwKgV2 Return molar masses [kg/mol] using PHREEQC BASIC GFW().

inSize = size(formulas);
formulas = string(formulas(:));
if isempty(formulas)
    error('phreeqcGfwKgV2:EmptyInput', 'Formula list must not be empty.');
end

mwKg = NaN(size(formulas));

dbPath = resolvePhreeqcDatabasePathV2(dbPath);
good = strlength(strtrim(formulas)) > 0;
if ~all(good)
    error('phreeqcGfwKgV2:BlankFormula', 'Blank formula passed to PHREEQC GFW().');
end

persistent gfwCache
if isempty(gfwCache)
    gfwCache = containers.Map('KeyType', 'char', 'ValueType', 'double');
end

keys = strings(size(formulas));
for i = 1:numel(formulas)
    keys(i) = dbPath + "|" + strtrim(formulas(i));
end

pendingMask = good;
for i = find(good(:)).'
    key = char(keys(i));
    if isKey(gfwCache, key)
        mwKg(i) = gfwCache(key);
        pendingMask(i) = false;
    end
end

if any(pendingMask)
    pending = formulas(pendingMask);
    uniquePending = unique(strtrim(pending), 'stable');
    valuesKg = eval_gfw_batch(uniquePending, dbPath);
    for i = 1:numel(uniquePending)
        key = char(dbPath + "|" + uniquePending(i));
        gfwCache(key) = valuesKg(i);
    end
    for i = find(pendingMask(:)).'
        key = char(keys(i));
        mwKg(i) = gfwCache(key);
    end
end

mwKg = reshape(mwKg, inSize);
end

function valuesKg = eval_gfw_batch(formulas, dbPath)
iph = get_iph_handle(dbPath);
formulas = string(formulas(:));
n = numel(formulas);

headings = "gfw_" + string(1:n);
expr = strings(1, n);
for i = 1:n
    expr(i) = "GFW(""" + replace(string(strtrim(formulas(i))), """", """""") + """)";
end

lines = strings(0, 1);
lines(end + 1, 1) = "DATABASE " + dbPath;
lines(end + 1, 1) = "";
lines(end + 1, 1) = "SOLUTION 1";
lines(end + 1, 1) = "    temp 25";
lines(end + 1, 1) = "    pH 7 charge";
lines(end + 1, 1) = "    units mol/kgw";
lines(end + 1, 1) = "    Na 1e-20";
lines(end + 1, 1) = "";
lines(end + 1, 1) = "SELECTED_OUTPUT 999";
lines(end + 1, 1) = "    -reset false";
lines(end + 1, 1) = "    -high_precision true";
lines(end + 1, 1) = "    -user_punch true";
lines(end + 1, 1) = "USER_PUNCH 999";
lines(end + 1, 1) = "    -headings " + strjoin(cellstr(headings.'), ' ');
lines(end + 1, 1) = "    -start";
lines(end + 1, 1) = "10 PUNCH " + strjoin(cellstr(expr), ', ');
lines(end + 1, 1) = "    -end";
lines(end + 1, 1) = "";
lines(end + 1, 1) = "RUN_CELLS";
lines(end + 1, 1) = "    -cells 1";
lines(end + 1, 1) = "END";
lines(end + 1, 1) = "";

err = iph.RunString(char(strjoin(lines, newline)));
if err ~= 0
    error('phreeqcGfwKgV2:RunFailed', 'PHREEQC GFW() failed: %s', iph.GetErrorString());
end

raw = iph.GetSelectedOutputValues(999);
arr = to_numeric_matrix(raw);
if isempty(arr)
    error('phreeqcGfwKgV2:EmptyOutput', 'PHREEQC GFW() returned no data.');
end

values = arr(end, :);
valuesKg = 1e-3 * values(:);
bad = ~(valuesKg > 0);
if any(bad)
    error('phreeqcGfwKgV2:InvalidMass', ...
        'PHREEQC GFW() returned invalid masses for: %s', strjoin(cellstr(formulas(bad)), ', '));
end
end

function iph = get_iph_handle(dbPath)
persistent iphCache
if isempty(iphCache)
    iphCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

key = char(dbPath);
if isKey(iphCache, key)
    iph = iphCache(key);
    return
end

iph = createIPhreeqcHandleV2(key);

iphCache(key) = iph;
end

function arr = to_numeric_matrix(values)
if isnumeric(values)
    arr = values;
    return
end

[nr, nc] = size(values);
arr = zeros(nr, nc);
for r = 1:nr
    for c = 1:nc
        v = values{r, c};
        if isnumeric(v)
            x = v;
        elseif isstring(v) || ischar(v)
            x = str2double(string(v));
        else
            x = v;
        end
        arr(r, c) = x;
    end
end
end
