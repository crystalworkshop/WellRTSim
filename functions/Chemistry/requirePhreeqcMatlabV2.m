function requirePhreeqcMatlabV2()
% requirePhreeqcMatlabV2 Ensure external PhreeqcMatlab is on the MATLAB path.

if exist('IPhreeqc', 'class') == 8 || ~isempty(which('IPhreeqc'))
    return
end

error('requirePhreeqcMatlabV2:MissingDependency', ...
    ['PhreeqcMatlab is not on the MATLAB path. Install it from ' ...
    'https://github.com/simulkade/PhreeqcMatlab and run its startup.m ' ...
    'before running WellRTSim.']);
end
