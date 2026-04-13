function state = ensureIPhreeqcHandleV2(state)
% ensureIPhreeqcHandleV2 Lazily create IPhreeqc handle if needed.

if ~isempty(state.chem.iph)
    return
end

dbPath = string(state.chem.db);
if strlength(dbPath) == 0
    dbPath = string(state.db);
end
dbPath = resolvePhreeqcDatabasePathV2(dbPath);

iph = createIPhreeqcHandleV2(dbPath);

state.chem.iph = iph;
state.chem.db = char(dbPath);
state.iphreeqc = 1;
end
