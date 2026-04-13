function setupPath = getChemistrySetupPath(state)
setupFile = "chemistry.md";

if isfield(state, 'chemistry_setup') && strlength(string(state.chemistry_setup)) > 0
    setupFile = string(state.chemistry_setup);
elseif isfield(state, 'chemistry_setup_file') && strlength(string(state.chemistry_setup_file)) > 0
    setupFile = string(state.chemistry_setup_file);
elseif isfield(state, 'chemistry_file') && strlength(string(state.chemistry_file)) > 0
    setupFile = string(state.chemistry_file);
end

setupFile = strtrim(setupFile);
if startsWith(setupFile, "/") || ~isempty(regexp(char(setupFile), '^[A-Za-z]:[\\/]', 'once'))
    setupPath = char(setupFile);
else
    setupPath = fullfile(state.SimDir, char(setupFile));
end
end
