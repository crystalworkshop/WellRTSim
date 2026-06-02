function state = initializeGraphics(state)
% Initialize all graphics components for the wellbore simulator
%
% Input:
%   state - State structure with basic parameters
%
% Output:
%   state - State structure with graphics components initialized

% Headless mode: build no UI at all. Install plain-struct stubs for the
% handles the main loop sets unconditionally (statusLabel + control buttons)
% and the control flags, so the simulation runs with zero graphics. Use this
% to diagnose render crashes that originate in the graphics layer, not the math.
if isfield(state, 'enablePlots') && ~state.enablePlots
    btnStub = struct('Enable', 'off');
    state.statusLabel = struct('Text', '');
    state.runBtn = btnStub;
    state.pauseBtn = btnStub;
    state.cancelBtn = btnStub;
    state.runFlag = true;
    state.pauseFlag = false;
    state.cancelFlag = false;
    fprintf('Graphics disabled (headless mode) - no UI created.\n');
    return;
end

% Suppress common UI component warnings
warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
warning('off', 'MATLAB:ui:javacomponent:FunctionToBeRemoved');

% MATLAB 2025a compatibility settings for VS Code
try
    % Force compatibility mode for UI figures
    if ~usejava('desktop')
        % In non-desktop mode, use compatible settings
        set(groot, 'defaultUIFigureVisible', 'on');
        set(groot, 'defaultUIAxesVisible', 'on');
    end
catch
    % Ignore if settings don't work in this MATLAB version
end

% Get screen size for proper window sizing
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% UI dimensions (0.8 × 0.8 of screen size, centered)
figWidth = round(0.8 * screenWidth);
figHeight = round(0.7 * screenHeight);

% Center the window on screen
figPosX = round((screenWidth - figWidth) / 2);
figPosY = round((screenHeight - figHeight) / 2);

tabMargin = 5;  % Reduced margin for smaller window

% Create a non-resizable main UI figure with fixed size
state.mainFig = uifigure('Position', [figPosX, figPosY, figWidth, figHeight], ...
                         'Name', 'Wellbore Simulator V2.0', ...
                         'Resize', 'on', ...
                         'WindowStyle', 'normal', ...
                         'AutoResizeChildren', 'on');

% Dark scheme colors
figBg   = [0.12 0.12 0.14];
panelBg = [0.16 0.16 0.19];
btnBg   = [0.22 0.22 0.26];
txtCol  = [0.94 0.94 0.96];
try
    state.mainFig.Color = figBg;
catch
end

% Root layout: plots above, controls below
state.mainLayout = uigridlayout(state.mainFig, [2,1]);
% Control bar height sized for 48px icons comfortably
state.mainLayout.RowHeight = {'1x', 48};
state.mainLayout.ColumnWidth = {'1x'};
state.mainLayout.Padding = [tabMargin tabMargin tabMargin tabMargin];

% Tab group in the top row
state.tabgp = uitabgroup(state.mainLayout);

% Create tabs
state.schematicTab = uitab(state.tabgp, 'Title', 'Well geometry');
state.schematicTab.AutoResizeChildren = 'off';
state.InitTab = uitab(state.tabgp, 'Title', 'Initial conditions');
state.TranTab = uitab(state.tabgp, 'Title', 'Profiles');
state.WellHeadTab = uitab(state.tabgp, 'Title', 'Well head');

% Controls layout in the bottom row
state.ctrlLayout = uigridlayout(state.mainLayout, [1,9]);
state.ctrlLayout.RowHeight = { '1x' };
state.ctrlLayout.ColumnWidth = {'1x','fit','fit','fit','fit','fit','fit','fit','fit'};
state.ctrlLayout.Padding = [8 8 8 8];
state.ctrlLayout.BackgroundColor = panelBg;

% Load and scale icons from Icons/ (2x)
try
    targetSize = [48 48];
    runIcon   = imresize(imread(fullfile('Icons','Play.png')),  targetSize, 'bicubic');
    pauseIcon = imresize(imread(fullfile('Icons','Pause.png')), targetSize, 'bicubic');
    stopIcon  = imresize(imread(fullfile('Icons','Stop.png')),  targetSize, 'bicubic');
    pdfIcon   = imresize(imread(fullfile('Icons','PDF.png')),   targetSize, 'bicubic');
    clearIcon = imresize(imread(fullfile('Icons','clear.png')), targetSize, 'bicubic');
    saveIcon  = imresize(imread(fullfile('Icons','1_139.png')), targetSize, 'bicubic');
catch
    % Fallback to paths if reading fails
    runIcon   = fullfile('Icons','Play.png');
    pauseIcon = fullfile('Icons','Pause.png');
    stopIcon  = fullfile('Icons','Stop.png');
    pdfIcon   = fullfile('Icons','PDF.png');
    clearIcon = fullfile('Icons','clear.png');
    saveIcon  = fullfile('Icons','1_139.png');
end

% Left: status label spans available space
state.statusLabel = uilabel(state.ctrlLayout, 'Text', 'Ready to start simulation...', ...
    'FontSize', 9, 'FontColor', txtCol, 'BackgroundColor', 'none');
state.statusLabel.Layout.Column = 1;

% Buttons in the center (icon-only)
state.runBtn = uibutton(state.ctrlLayout, 'push', 'Text','', 'Icon', runIcon, ...
    'Tooltip','Run', 'BackgroundColor', btnBg);
state.runBtn.Layout.Column = 2;

% Plot frequency spinner (next to Run)
state.plotFreqSpinner = uispinner(state.ctrlLayout, ...
    'Limits', [1 Inf], 'Value', 1, 'Step', 1, ...
    'Tooltip', 'Plot every N steps', ...
    'BackgroundColor', btnBg);
state.plotFreqSpinner.Layout.Column = 3;

state.pauseBtn = uibutton(state.ctrlLayout, 'push', 'Text','', 'Icon', pauseIcon, ...
    'Tooltip','Pause', 'BackgroundColor', btnBg);
state.pauseBtn.Layout.Column = 4;

state.savePhreeqcBtn = uibutton(state.ctrlLayout, 'push', 'Text','PHR', ...
    'Tooltip','Save latest PHREEQC input (.phr)', 'BackgroundColor', btnBg);
state.savePhreeqcBtn.Layout.Column = 5;

state.cancelBtn = uibutton(state.ctrlLayout, 'push', 'Text','', 'Icon', stopIcon, ...
    'Tooltip','Stop', 'BackgroundColor', btnBg);
state.cancelBtn.Layout.Column = 6;

% Snapshot button (export current UI to PDF)
state.saveBtn = uibutton(state.ctrlLayout, 'push', 'Text','', 'Icon', pdfIcon, ...
    'Tooltip','Snapshot (PDF)', 'BackgroundColor', btnBg);
state.saveBtn.Layout.Column = 7;

% Save results (MAT file)
state.saveResBtn = uibutton(state.ctrlLayout, 'push', 'Text','', 'Icon', saveIcon, ...
    'Tooltip','Save results (MAT)', 'BackgroundColor', btnBg);
state.saveResBtn.Layout.Column = 8;

% Clear plots button
state.clearBtn = uibutton(state.ctrlLayout, 'push', 'Text','', 'Icon', clearIcon, ...
    'Tooltip','Clear plots', 'BackgroundColor', btnBg);
state.clearBtn.Layout.Column = 9;

% Axes for transient and wellhead plots
state.TranAxes = setupTransientAxes(state.TranTab, state.pressureUnitLabel);
state.WellAxes = setupWellheadAxes(state.WellHeadTab, state.pressureUnitLabel, ...
    state.chem.plotLabels, state.chem.gasMask);

% Restore the dedicated well-geometry view on its own tab.
wellName = "Well";
if isfield(state, 'wellData') && ~isempty(state.wellData) && ...
        ismember('WellName', state.wellData.Properties.VariableNames) && height(state.wellData) >= 1
    wellName = string(state.wellData.WellName(1));
end
createWellPlots(state, state.wellData, state.deviationData, 0, state.Lp, char(wellName));

state = configureControlCallbacks(state);

% Force MATLAB to update the UI layout
drawnow;

fprintf('Graphics components initialized successfully.\n');
end

function state = configureControlCallbacks(state)
state.runBtn.ButtonPushedFcn = @runCallback;
state.pauseBtn.ButtonPushedFcn = @pauseCallback;
state.cancelBtn.ButtonPushedFcn = @cancelCallback;
state.savePhreeqcBtn.ButtonPushedFcn = @savePhreeqcCallback;
state.plotFreqSpinner.ValueChangedFcn = @plotFreqChanged;
state.plotFreqSpinner.ValueChangingFcn = @plotFreqChanging;
state.saveBtn.ButtonPushedFcn = @saveSnapshotCallback;
state.clearBtn.ButtonPushedFcn = @clearPlotsCallback;
state.saveResBtn.ButtonPushedFcn = @saveResultsCallback;

state.runBtn.Enable = 'off';
state.pauseBtn.Enable = 'on';
state.cancelBtn.Enable = 'on';
state.pauseFlag = false;
state.cancelFlag = false;
state.runFlag = true;
state.plotFreqSpinner.Value = max(1, round(state.pltf));

userData = struct('TranTab', state.TranTab);
state.runBtn.UserData = userData;
state.pauseBtn.UserData = userData;
state.cancelBtn.UserData = userData;
end

function runCallback(~, ~)
state = getBaseState();
state.pauseFlag = false;
state.runFlag = true;
setBaseState(state);
end

function pauseCallback(~, ~)
state = getBaseState();
state.pauseFlag = true;
setBaseState(state);
end

function cancelCallback(src, ~)
state = getBaseState();
state.cancelFlag = true;
state.runFlag = false;
setBaseState(state);

try
    ud = src.UserData;
    if isstruct(ud)
        refreshTransientPlots(state, state.Y0, 'Stopped: showing last Y0');
        tranTab = ud.TranTab;
        tranTab.Parent.SelectedTab = tranTab;
    end
catch
    % Fall back to the main-loop stop handler if direct plotting fails.
end
end

function plotFreqChanged(src, ~)
state = getBaseState();
state.pltf = max(1, round(src.Value));
state.statusLabel.Text = sprintf('Plot every %d steps', state.pltf);
setBaseState(state);
end

function plotFreqChanging(~, evt)
state = getBaseState();
state.pltf = max(1, round(evt.Value));
state.statusLabel.Text = sprintf('Plot every %d steps', state.pltf);
setBaseState(state);
end

function savePhreeqcCallback(~, ~)
state = getBaseState();
try
    scriptText = string(state.chem.lastPhreeqcScript);
    if strlength(scriptText) == 0
        state.statusLabel.Text = 'No PHREEQC script available yet';
        setBaseState(state);
        return;
    end

    ensureDirectory(state.SimDir);
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
    outFile = fullfile(state.SimDir, ['phreeqc_saved_' stamp '.phr']);
    writeTextFile(outFile, char(scriptText));
    state.statusLabel.Text = ['PHREEQC file saved: ' outFile];
    setBaseState(state);
catch ME
    warning('initializeGraphics:SavePhreeqcFailed', '%s', ME.message);
end
end

function saveSnapshotCallback(~, ~)
state = getBaseState();
try
    state.saveCounter = state.saveCounter + 1;
    ensureDirectory(state.SimDir);
    outFile = fullfile(state.SimDir, sprintf('result_%03d.pdf', state.saveCounter));
    exportTransientToPDF(state, outFile);
    state.statusLabel.Text = ['Snapshot saved: ' outFile];
    setBaseState(state);
catch ME
    warning('initializeGraphics:SaveSnapshotFailed', '%s', ME.message);
end
end

function clearPlotsCallback(~, ~)
state = getBaseState();
try
    fns = fieldnames(state.TranAxes);
    for k = 1:numel(fns)
        ax = state.TranAxes.(fns{k});
        if isvalid(ax)
            delete(findall(ax, 'Type', 'line'));
            delete(findall(ax, 'Type', 'text', 'Tag', 'HoverDataTipText'));
        end
    end
    state.statusLabel.Text = 'Plots cleared';
    setBaseState(state);
catch ME
    warning('initializeGraphics:ClearPlotsFailed', '%s', ME.message);
end
end

function saveResultsCallback(~, ~)
state = getBaseState();
try
    state.resultsCounter = state.resultsCounter + 1;
    ensureDirectory(state.SimDir);
    outFile = fullfile(state.SimDir, sprintf('results_%03d.mat', state.resultsCounter));
    t_history = evalin('base', 't_history');
    Y_history = evalin('base', 'Y_history');
    try
        save(outFile, 'state', 't_history', 'Y_history', '-v7.3');
    catch
        save(outFile, 'state', 't_history', 'Y_history');
    end
    state.statusLabel.Text = ['Saved results: ' outFile];
    setBaseState(state);
catch ME
    warning('initializeGraphics:SaveResultsFailed', '%s', ME.message);
end
end

function state = getBaseState()
state = evalin('base', 'state');
end

function setBaseState(state)
assignin('base', 'state', state);
end

function ensureDirectory(pathStr)
if exist(pathStr, 'dir') ~= 7
    mkdir(pathStr);
end
end

function writeTextFile(path, content)
fid = fopen(path, 'w');
if fid < 0
    error('initializeGraphics:CannotWriteFile', ...
        'Cannot open file for writing: %s', path);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s', content);
clear cleaner;
end
