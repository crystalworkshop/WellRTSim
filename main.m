% WellRTSim - Wellbore Reactive Transport Simulator
% Copyright (c) 2026  Oleg Melnik
% Licensed under CC BY-NC-SA 4.0
% For commercial use, contact: oleg,melnik@earth.ox.ac.uk, oemelnik@gmail.com
%%
close all force
clearvars

% Initialize state structure with parameters
addpath('functions/Chemistry');
addpath('functions/Hydrodynamics');
addpath('functions/Wellspec');
addpath('functions/Graphics');
addpath('functions/IO');
addpath('functions/init');

state = struct(); % Start with empty state
SimDir='Simulation/Krafla/KJ-9/';
% SimDir='Simulation/Tonkin/LO_T_heat/';
SimFile='params.md';
%SimDir='Simulation/WellSim/RK-5_id1/';
%SimDir='Simulation/Dixie_Valley/well_84-7/';
state.SimDir = SimDir;
state.SimFile = SimFile;

% Set false to disable all live plotting/drawnow (diagnose graphics-render crashes)
state.enablePlots = true;

state = initializeState(state);
state = initializeGraphics(state);

state = initializeSteadyState(state); 
if state.iBC_top == 3
    state = updateTopPressureFromWHP(state);
end
% Initialize chemistry initial conditions using current (P,T) profile

% Print simulation parameters
fprintf('Running flow simulation:\n');
% Initialize HDF5 results file in the simulation directory
state = initResultsH5(state);

if state.calc_chem==1
    state = chemistryInitializeStateV2(state);
    state.wellhead.CoutPpm = liquidMassfracToPpmV2(state.C(:, end));
end
if state.enablePlots
    updateWellheadPlots(state);
    refreshTransientPlots(state, state.Y, 'Initial conditions');
    drawnow;
end
[k,n] = size(state.Y);
% For storing results at different time steps
num_saves = 10;  % Number of time steps to save
save_interval = ceil(state.tfin / (num_saves - 1));
Y_history = zeros(k, n, num_saves);
t_history = zeros(1, num_saves);
Y_history(:,:,1) = state.Y;
t_history(1) = 0;
save_idx = 2;

% Time stepping parameters (seconds)

fprintf('Starting simulation...\n');
if state.calc_chem == 1 && state.stat_chem > 0
    fprintf('Chemical transport starts at t = %.3f s.\n', state.stat_chem);
end
fprintf('Time = %7.3f, dt = %7.3f\n', state.tt, state.dt);
state.stepIndex = 2;

%% Time loop
% Switch to Transient tab and plot initial state once (keep previous)
if state.enablePlots
    refreshTransientPlots(state, state.Y, ['Time = ' num2str(state.tt/3600/24) ' days']);
    state.tabgp.SelectedTab = state.TranTab;
end
state = appendProfileH5(state, true);
istep = 0; % time step counter

while state.tt <= state.tfin && state.runFlag
    % Allow UI callbacks to execute during the loop
    if state.enablePlots, drawnow limitrate; end
    % Check for pause
    if state.pauseFlag
        state.statusLabel.Text = 'Simulation paused - showing Y0 values';
        % Plot Y0 values in transient tab when paused (overlay)
        if state.enablePlots
            refreshTransientPlots(state, state.Y0, ['Paused @ ' num2str(state.tt/3600/24) ' days']);
            state.tabgp.SelectedTab = state.TranTab;
        end

        % Enable run button while paused
        state.runBtn.Enable = 'on';
        state.pauseBtn.Enable = 'off';

        % Wait for resume
        while state.pauseFlag && ~state.cancelFlag && state.runFlag
            pause(0.1);
            drawnow;
        end

        if state.cancelFlag, break; end
        if ~state.runFlag, break; end

        % Re-enable pause, disable run
        state.pauseBtn.Enable = 'on';
        state.runBtn.Enable = 'off';
        state.statusLabel.Text = sprintf('Resuming simulation... Time %f s', state.tt);
    end

    % Check if cancelled
    if state.cancelFlag, break; end

    % Take one time step
    state = updateTopPressureFromWHP(state);
    state.Dpprev=state.Dp;
    [Y, tol1, tol2,it, state] = OneStep(state);
    istep = istep + 1;

    if state.cancelFlag
        break;
    end
    if ~state.stepConverged
        fprintf('Simulation did not converge at time t = %7.3f with dt = %7.3e\n', state.tt, state.dt);
        state.cancelFlag = true;
        break;
    end

    state.wellhead.iterations = it;
    state.wellhead.tol1 = tol1;
    
    % state.P_top=max(2e6,state.P_top-1e5*state.tt);

    % Update status display
    progress = min(1, state.tt/state.tfin);
    state.statusLabel.Text = sprintf('Running: Time %.3f s (%.1f%% complete)', state.tt, progress*100);
    if state.enablePlots, drawnow limitrate; end % process UI events and keep app responsive

    if state.cancelFlag, break; end
    state = updateWellheadDiagnostics(state);
    %% chemistry model
    chemTimeS = state.tt + state.dt;
    if state.calc_chem == 1 && chemTimeS >= state.stat_chem
        state = chemistryStepV2(state);
        state.chemSample.iterations = 1;
        state.chemSample.tol1 = 0;
        state = appendChemistryH5(state);
    end
    if state.calc_chem == 1
        state.wellhead.CoutPpm = liquidMassfracToPpmV2(state.C(:, end));
    end
    state = appendWellheadH5(state); % Save wellhead time-series after outlet chemistry is updated
    if state.enablePlots, updateWellheadPlots(state); end

    state.stepIndex = state.stepIndex + 1;
    state.Y0 = Y;
    state.Y = Y;
    state.tt = state.tt + state.dt;

    % Plot every 'pltf' steps and keep previous lines
    if state.enablePlots && mod(istep, max(1, round(state.pltf))) == 0
        refreshTransientPlots(state, state.Y, ['t = ' num2str(state.tt/3600/24) ' days']);
        updateWellheadPlots(state);
        state.tabgp.SelectedTab = state.TranTab;
        drawnow limitrate;
    end
    % if state.tt>3600*24*10, state.pltf=1; end
    % Save spatial profiles to HDF5/CSV at configured cadence
    state = appendProfileH5(state);

    % Output current time and error norms
    pScale = state.pressureUnitScale;
    pLabel = char(state.pressureUnitLabel);
    fprintf('Time = %7.3f, dt = %7.3e, iters = %d, tol1 = %10.4e, tol2 = %10.4e P_bot = %10.4e %s\n', ...
        state.tt, state.dt, it, tol1, tol2, state.Y(1,1)/pScale, pLabel);

    % Save history at specified intervals
    if state.tt >= (save_idx-1)*save_interval && save_idx <= num_saves
        Y_history(:,:,save_idx) = state.Y;
        t_history(save_idx) = state.tt;
        save_idx = save_idx + 1;
    end

    % Check for convergence
    if ~isfinite(tol1)
        fprintf('Simulation did not converge at time t = %7.3f\n', state.tt);
        state.cancelFlag = true;
        break;
    end

    % Adopt the controller's next-step dt (OneStep left state.dt at the dt that
    % was actually integrated, so the time advance above used the solved step).
    state.dt = state.dtNext;
end

% Final results and cleanup
if state.cancelFlag
    state.statusLabel.Text = 'Simulation stopped';
    if state.enablePlots
        refreshTransientPlots(state, state.Y0, ['Final t = ' num2str(state.tt/3600/24) ' days']);
        state.tabgp.SelectedTab = state.TranTab;
    end
else
    state.statusLabel.Text = 'Simulation complete!';
    fprintf('\nSimulation completed. Plotting results...\n');
    if state.enablePlots
        refreshTransientPlots(state, state.Y, ['Final t = ' num2str(state.tt/3600/24) ' days']);
        state.tabgp.SelectedTab = state.TranTab;
    end
end

% Reset control buttons
state.runBtn.Enable = 'on';
state.pauseBtn.Enable = 'off';
state.cancelBtn.Enable = 'off';
