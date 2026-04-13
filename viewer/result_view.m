function result_view(targetPath)
% WellRTSim - Wellbore Reactive Transport Simulator
% Copyright (c) 2026  Oleg Melnik
% Licensed under CC BY-NC-SA 4.0
% For commercial use, contact: oleg,melnik@earth.ox.ac.uk, oemelnik@gmail.com
% result_view Launch the standalone results viewer for simulation outputs.
%
% Usage:
%   addpath(fullfile(pwd, 'viewer'));
%   result_view
%   result_view(simDir)
%   result_view(resultsFile)

    if nargin < 1
        targetPath = "";
    end

    toolDir = fileparts(mfilename('fullpath'));
    rv_launchResultView(targetPath, toolDir);
end
