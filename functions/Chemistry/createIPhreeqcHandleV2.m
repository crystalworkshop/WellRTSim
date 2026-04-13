function iph = createIPhreeqcHandleV2(dbPath)
% createIPhreeqcHandleV2 Create a configured IPhreeqc handle.

dbPath = char(string(dbPath));

try
    iph = IPhreeqc();
catch ME
    error('createIPhreeqcHandleV2:IPhreeqcMissing', ...
        'Unable to create IPhreeqc instance: %s', ME.message);
end

status = iph.LoadDatabase(dbPath);
if status ~= 0
    msg = strtrim(string(iph.GetErrorString()));
    if strlength(msg) == 0
        msg = "unknown error";
    end
    error('createIPhreeqcHandleV2:DatabaseLoadFailed', ...
        'PHREEQC failed to load database "%s": %s', dbPath, char(msg));
end

safe_iph_call(iph, 'SetOutputStringOn', true);
safe_iph_call(iph, 'SetErrorStringOn', true);
safe_iph_call(iph, 'SetSelectedOutputFileOn', false);
safe_iph_call(iph, 'SetOutputFileOn', false);
safe_iph_call(iph, 'SetLogFileOn', false);
safe_iph_call(iph, 'SetDumpFileOn', false);
end

function safe_iph_call(obj, methodName, varargin)
if ismethod(obj, methodName)
    obj.(methodName)(varargin{:});
end
end
