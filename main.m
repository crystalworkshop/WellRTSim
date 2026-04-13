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
% SimDir='Simulation/Sumatra/AA03_1380m/';
SimDir='Simulation/Krafla/KJ-9/';
SimFile='params.md';
%SimDir='Simulation/WellSim/RK-5_id1/';
%SimDir='Simulation/Dixie_Valley/well_84-7/';
state.SimDir = SimDir;
state.SimFile = SimFile;
state = initializeState(state);
state = initializeGraphics(state);

state = initializeSteadyState(state); 
% Initialize chemistry initial conditions using current (P,T) profile

% Print simulation parameters
fprintf('Running flow simulation:\n');
% Initialize HDF5 results file in the simulation directory
state = initResultsH5(state);

if state.calc_chem==1
    state = chemistryInitializeStateV2(state);
    state.CaInKg(1) = 0;
    state.CaOutKg(1) = 0;
    state.CaPrecipKg(1) = 0;
    state.CaFluidDeltaKg(1) = 0;
    state.CaBalanceResidualKg(1) = 0;
end
state.CoutPpm(:, 1) = liquidMassfracToPpmV2(state.C(:, end));
updateWellheadPlots(state);

refreshTransientPlots(state, state.Y, 'Initial conditions');
drawnow;
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
refreshTransientPlots(state, state.Y, ['Time = ' num2str(state.tt/3600/24) ' days']);
state.tabgp.SelectedTab = state.TranTab;
state = appendProfileH5(state, true);
istep = 0; % time step counter

while state.tt <= state.tfin && state.runFlag
    % Allow UI callbacks to execute during the loop
    drawnow limitrate;
    % Check for pause
    if state.pauseFlag
        state.statusLabel.Text = 'Simulation paused - showing Y0 values';
        % Plot Y0 values in transient tab when paused (overlay)
        refreshTransientPlots(state, state.Y0, ['Paused @ ' num2str(state.tt/3600/24) ' days']);
        state.tabgp.SelectedTab = state.TranTab;

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
    state.Dpprev=state.Dp;
    [Y, tol1, tol2,it, state] = OneStep(state);
    istep = istep + 1;
    state = updateTopPressureFromWHP(state);
    idx = state.stepIndex;
    state.Iterations(idx) = it;
    state.Tol1(idx) = tol1;
    
    % state.P_top=max(2e6,state.P_top-1e5*state.tt);

    % Update status display
    progress = min(1, state.tt/state.tfin);
    state.statusLabel.Text = sprintf('Running: Time %.3f s (%.1f%% complete)', state.tt, progress*100);
    drawnow limitrate; % process UI events and keep app responsive

    if state.cancelFlag, break; end
    state = updateWellheadDiagnostics(state, idx, state.Y(3, end - 1), state.Y(3, end));
    %% chemistry model
    chemTimeS = state.tt + state.dt;
    if state.calc_chem == 1 && chemTimeS >= state.stat_chem
        state = chemistryStepV2(state);
        state.ChemIterations(idx) = 1;
        state.ChemTol1(idx) = 0;
        caBal = state.chem.lastCalciumBalance;
        state.CaInKg(idx) = caBal.inKg;
        state.CaOutKg(idx) = caBal.outKg;
        state.CaPrecipKg(idx) = caBal.precipHydroKg;
        state.CaFluidDeltaKg(idx) = caBal.fluidDeltaKg;
        state.CaBalanceResidualKg(idx) = caBal.residualKg;
        fprintf('  Ca hydro balance [kg/dt]: in=%10.4e, precip=%10.4e, out=%10.4e, dfluid=%10.4e, resid=%10.4e\n', ...
            caBal.inKg, caBal.precipHydroKg, caBal.outKg, caBal.fluidDeltaKg, caBal.residualKg);
        fprintf('  Ca chem  balance [kg/flight]: pre=%10.4e, precip=%10.4e, post=%10.4e, resid=%10.4e\n', ...
            caBal.chemInKg, caBal.precipKg, caBal.chemOutKg, caBal.chemResidualKg);
        state = appendChemistryH5(state);
    end
    state.CoutPpm(:, idx) = liquidMassfracToPpmV2(state.C(:, end));
    state = appendWellheadH5(state); % Save wellhead time-series after outlet chemistry is updated
    updateWellheadPlots(state);

    state.stepIndex = state.stepIndex + 1;
    state.Y0 = Y;
    state.Y = Y;
    state.tt = state.tt + state.dt;

    % Plot every 'pltf' steps and keep previous lines
    if mod(istep, max(1, round(state.pltf))) == 0
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
    if isnan(tol1)
        fprintf('Simulation did not converged at time t = %7.3f\n', state.tt);
        state.cancelFlag = true;
        break;
    end
end

% Final results and cleanup
if state.cancelFlag
    state.statusLabel.Text = 'Simulation stopped';
    refreshTransientPlots(state, state.Y0, ['Final t = ' num2str(state.tt/3600/24) ' days']);
    state.tabgp.SelectedTab = state.TranTab;
else
    state.statusLabel.Text = 'Simulation complete!';
    fprintf('\nSimulation completed. Plotting results...\n');
    refreshTransientPlots(state, state.Y, ['Final t = ' num2str(state.tt/3600/24) ' days']);
    state.tabgp.SelectedTab = state.TranTab;
end

% Reset control buttons
state.runBtn.Enable = 'on';
state.pauseBtn.Enable = 'off';
state.cancelBtn.Enable = 'off';
