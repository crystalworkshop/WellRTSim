function cfg = parseChemistrySetup(mdPath)
%PARSECHEMISTRYSETUP Parse chemistry markdown setup file.
%   cfg = parseChemistrySetup(mdPath) parses chemistry.md and returns a
%   normalized struct used by PHREEQC input generation.

if nargin < 1 || strlength(string(mdPath)) == 0
    error('parseChemistrySetup requires a path to chemistry markdown file.');
end

[mdPath, simDir] = resolve_setup_file(mdPath);
secs = read_md_sections(mdPath);

aq = parse_aqueous(get_section_lines(secs, "AQUEOUS SPECIES"));
gss = parse_gases(get_section_lines(secs, "GASES"));
mins = parse_minerals(get_section_lines(secs, "MINERALS"));
tpl = parse_template_file(get_section_lines(secs, "TEMPLATE FILE"));
db = parse_database_path(get_section_lines(secs, "DATABASE FILE PATH"), mdPath);

if strlength(tpl.file_name) == 0
    tpl.file_name = "chemistry.pht";
end

tplPath = string(tpl.file_name);
if ~is_absolute_path(tplPath)
    tplPath = fullfile(simDir, tplPath);
end
tplPath = string(char(java.io.File(char(tplPath)).getCanonicalPath()));

cfg = struct();
cfg.paths = struct();
cfg.paths.simulation_dir = string(simDir);
cfg.paths.setup_md = string(mdPath);
cfg.paths.template_path = tplPath;

cfg.template_file = string(tpl.file_name);
cfg.template_path = tplPath;
cfg.generate_phreeqc = tpl.generate_phreeqc;
cfg.database_path = string(db.path);

cfg.aqueous_names = aq.names;
cfg.aqueous_suffix = aq.suffix;
cfg.aqueous_concentration_ppm = aq.concentration_ppm;
cfg.aqueous_concentration_molal = aq.concentration_molal;
cfg.gas_names = gss.names;
cfg.gas_partition = gss.partition;

cfg.kinetics_names = mins.kinetics;
cfg.mineral_formula = mins.formula;
cfg.mineral_tau_s = mins.tau_s;
cfg.mineral_density_kg_m3 = mins.density_kg_m3;
end

%% ========================= Path Helpers ===========================
function [setupFile, simDir] = resolve_setup_file(inputPath)
if isstring(inputPath)
    inputPath = char(inputPath);
end
if ~ischar(inputPath)
    error('chemistry setup path must be char or string.');
end

raw = strtrim(inputPath);
if isempty(raw)
    error('chemistry setup path is empty.');
end

if exist(raw, 'dir') == 7
    setupFile = fullfile(raw, 'chemistry.md');
elseif exist(raw, 'file') == 2
    setupFile = raw;
else
    error('Chemistry setup path not found: %s', raw);
end

if exist(setupFile, 'file') ~= 2
    error('Chemistry markdown file not found: %s', setupFile);
end

setupFile = char(java.io.File(setupFile).getCanonicalPath());
simDir = fileparts(setupFile);
end

%% ====================== Markdown Processing =======================
function secs = read_md_sections(setupFile)
txt = fileread(setupFile);
lines = splitlines(string(txt));

secs = struct('title', {}, 'lines', {});
iSec = 0;
for i = 1:numel(lines)
    line = string(lines(i));
    tok = regexp(char(line), '^\s*##\s*(.+?)\s*$', 'tokens', 'once');
    if ~isempty(tok)
        iSec = iSec + 1;
        secs(iSec).title = string(tok{1}); %#ok<AGROW>
        secs(iSec).lines = strings(0, 1); %#ok<AGROW>
        continue
    end
    if iSec == 0
        continue
    end
    secs(iSec).lines(end + 1, 1) = line; %#ok<AGROW>
end

if isempty(secs)
    error('No level-2 sections (## ...) found in %s', setupFile);
end
end

function lines = get_section_lines(secs, titleKey)
lines = strings(0, 1);
if isempty(secs)
    return
end

keyNeed = norm_key(titleKey);
for i = 1:numel(secs)
    if norm_key(secs(i).title) == keyNeed
        lines = secs(i).lines;
        return
    end
end
end

function out = clean_payload(lines)
out = strings(0, 1);
for i = 1:numel(lines)
    raw = char(lines(i));
    rawTrim = strtrim(raw);
    if isempty(rawTrim)
        continue
    end
    if startsWith(rawTrim, '#')
        continue
    end
    noTailComment = regexprep(raw, '\s+#.*$', '');
    line = strtrim(noTailComment);
    if isempty(line)
        continue
    end
    out(end + 1, 1) = string(line); %#ok<AGROW>
end
end

function k = norm_key(v)
k = lower(strtrim(string(v)));
k = regexprep(k, '[^a-z0-9]', '');
end

%% ======================= Section Parsers ==========================
function aq = parse_aqueous(lines)
payload = clean_payload(lines);

names = strings(0, 1);
suffixRaw = strings(0, 1);
ppmRaw = NaN(0, 1);
molalRaw = NaN(0, 1);
pendingSuffix = false;

for i = 1:numel(payload)
    line = string(payload(i));

    if pendingSuffix
        suffixRaw = parse_text_list(line);
        pendingSuffix = false;
        continue
    end

    [ok, key, val] = split_kv(line);
    if ok
        k = norm_key(key);
        if any(k == ["suffix", "suffixes", "aqsuffix", "aqueoussuffix"])
            suffixRaw = parse_text_list(val);
        elseif any(k == ["concentrationppm", "concppm", "ppm"])
            ppmRaw = parse_num_list(val);
        elseif any(k == ["concentrationmolal", "concentrationmolkgw", "molal", "molkgw"])
            molalRaw = parse_num_list(val);
        end
        continue
    end

    if any(norm_key(line) == ["suffix", "suffixes"])
        pendingSuffix = true;
        continue
    end

    if contains(line, ':')
        continue
    end

    names = [names; split_name_list(line)]; %#ok<AGROW>
end

if isempty(names)
    error('Aqueous species section is empty.');
end

n = numel(names);
suffix = pad_text_vector(suffixRaw, n, "");
ppm = pad_num_vector(ppmRaw, n, NaN);
molal = pad_num_vector(molalRaw, n, NaN);

aq = struct();
aq.names = names;
aq.suffix = suffix;
aq.concentration_ppm = ppm;
aq.concentration_molal = molal;
end

function mins = parse_minerals(lines)
payload = clean_payload(lines);

kin = strings(0, 1);
formula = strings(0, 1);
tau = NaN(0, 1);
density = NaN(0, 1);

for i = 1:numel(payload)
    [ok, key, val] = split_kv(payload(i));
    if ~ok
        continue
    end
    k = norm_key(key);
    if any(k == ["kinetics", "kinetic"])
        kin = split_name_list(val);
    elseif any(k == ["formula", "formulas"])
        formula = split_name_list(val);
    elseif any(k == ["taus", "tau", "tausseconds", "taus"])
        tau = parse_num_list(val);
    elseif any(k == ["density", "densitykgm3", "densitykgperm3"])
        density = parse_num_list(val);
    end
end

if isempty(kin)
    error('Minerals section: Kinetics list is required.');
end

nk = numel(kin);
formula = pad_text_vector(formula, nk, "");
tau = pad_num_vector(tau, nk, NaN);
density = pad_num_vector(density, nk, NaN);

mins = struct();
mins.kinetics = kin;
mins.formula = formula;
mins.tau_s = tau;
mins.density_kg_m3 = density;
end

function gss = parse_gases(lines)
payload = clean_payload(lines);
names = strings(0, 1);
partition = 2;

for i = 1:numel(payload)
    line = string(payload(i));
    [ok, key, val] = split_kv(line);
    if ok
        k = norm_key(key);
        if any(k == ["gases", "gas", "gasspecies", "gasphase"])
            names = [names; split_name_list(val)]; %#ok<AGROW>
        elseif any(k == ["partition", "gaspartition", "partitionmode"])
            partVal = str2double(string(val));
            if isfinite(partVal)
                partition = round(partVal);
            end
        end
        continue
    end
    if contains(line, ':')
        continue
    end
    names = [names; split_name_list(line)]; %#ok<AGROW>
end

if isempty(names)
    gss.names = strings(0, 1);
else
    names = strtrim(names);
    names = names(strlength(names) > 0);
    gss.names = unique(names, 'stable');
end

if ~ismember(partition, [1, 2])
    error('Gases section: partition must be 1 (PHREEQC gas phase) or 2 (MATLAB analytical model).');
end
gss.partition = partition;
end

function tpl = parse_template_file(lines)
payload = clean_payload(lines);
tpl = struct();
tpl.file_name = "";
tpl.generate_phreeqc = true;

for i = 1:numel(payload)
    line = string(payload(i));
    [ok, key, val] = split_kv(line);
    if ok
        if norm_key(key) == "generatephreeqc"
            t = lower(strtrim(string(val)));
            if any(t == ["1", "true", "yes", "y", "on"])
                tpl.generate_phreeqc = true;
            elseif any(t == ["0", "false", "no", "n", "off"])
                tpl.generate_phreeqc = false;
            else
                tpl.generate_phreeqc = true;
            end
        end
        continue
    end
    if strlength(tpl.file_name) == 0
        tpl.file_name = strip_wrapping_quotes(line);
    end
end
end

function db = parse_database_path(lines, setupFile)
payload = clean_payload(lines);
db = struct();
db.path = "";

for i = 1:numel(payload)
    line = strtrim(string(payload(i)));
    [ok, key, val] = split_kv(line);
    if ok
        k = norm_key(key);
        if any(k == ["database", "databasefile", "db", "path", "databasefilepath"])
            db.path = strip_wrapping_quotes(strtrim(string(val)));
            break
        end
    else
        db.path = strip_wrapping_quotes(line);
        break
    end
end

if strlength(db.path) == 0
    error('Database path is missing in %s.', setupFile);
end
end

%% ============================ Utils ===============================
function [ok, key, val] = split_kv(line)
if ~isempty(regexp(char(strtrim(string(line))), '^[A-Za-z]:[\\/]', 'once'))
    ok = false;
    key = "";
    val = "";
    return
end

tok = regexp(char(line), '^\s*([^:]+?)\s*:\s*(.*?)\s*$', 'tokens', 'once');
if isempty(tok)
    ok = false;
    key = "";
    val = "";
    return
end
ok = true;
key = string(tok{1});
val = string(tok{2});
end

function list = split_name_list(v)
parts = split(string(v), ',');
list = strings(0, 1);
for i = 1:numel(parts)
    chunk = strtrim(parts(i));
    if strlength(chunk) == 0
        continue
    end
    tokens = split(chunk);
    for j = 1:numel(tokens)
        nm = strip_wrapping_quotes(strtrim(tokens(j)));
        nm = regexprep(nm, '[\.;]+$', '');
        if strlength(nm) == 0
            continue
        end
        list(end + 1, 1) = nm; %#ok<AGROW>
    end
end
end

function list = parse_text_list(v)
txt = string(v);
if contains(txt, ',')
    parts = split(txt, ',');
    list = strings(numel(parts), 1);
    for i = 1:numel(parts)
        tok = strip_wrapping_quotes(strtrim(string(parts(i))));
        tok = regexprep(tok, '[\.;]+$', '');
        list(i) = tok;
    end
else
    parts = split(txt);
    list = strings(0, 1);
    for i = 1:numel(parts)
        tok = strip_wrapping_quotes(strtrim(string(parts(i))));
        tok = regexprep(tok, '[\.;]+$', '');
        if strlength(tok) == 0
            continue
        end
        list(end + 1, 1) = tok; %#ok<AGROW>
    end
end
end

function vals = parse_num_list(v)
txt = char(string(v));
tok = regexp(txt, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
if isempty(tok)
    vals = NaN(0, 1);
    return
end
vals = str2double(string(tok(:)));
end

function out = pad_text_vector(in, n, fillValue)
if nargin < 3
    fillValue = "";
end
out = strings(max(0, n), 1);
if isempty(in)
    out(:) = string(fillValue);
    return
end
in = string(in(:));
nCopy = min(n, numel(in));
if nCopy > 0
    out(1:nCopy) = in(1:nCopy);
end
if nCopy < n
    out(nCopy+1:n) = string(fillValue);
end
end

function out = pad_num_vector(in, n, fillValue)
if nargin < 3
    fillValue = NaN;
end
out = fillValue * ones(max(0, n), 1);
if isempty(in)
    return
end
in = in(:);
nCopy = min(n, numel(in));
if nCopy > 0
    out(1:nCopy) = in(1:nCopy);
end
end

function out = strip_wrapping_quotes(in)
out = string(in);
out = regexprep(out, '^[\"'']+|[\"'']+$', '');
end

function tf = is_absolute_path(p)
s = char(strtrim(string(p)));
if isempty(s)
    tf = false;
    return
end
if ispc
    tf = ~isempty(regexp(s, '^[A-Za-z]:[\\/]|^\\\\', 'once'));
else
    tf = startsWith(string(s), "/");
end
end
