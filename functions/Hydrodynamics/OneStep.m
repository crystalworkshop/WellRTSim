function [Y, tol1, tol2, it, state] = OneStep(state)
% Advance one timestep with Newton-Raphson. If the step fails to converge
% within maxiter iterations, reject it and retry with a smaller dt (down to
% dtMin). The next-step dt is set from the iteration count:
%   it < 10  -> grow   (dt * dt_increment)
%   it > 20  -> shrink (dt / dt_increment)
%   otherwise -> hold.
% state.dt is left at the dt that was actually integrated this step, so the
% caller advances time by the solved step; state.dtNext carries the next dt.
toler = state.epsQ;
Y0 = state.Y0;
dtMin = 1e-6;
state.stepConverged = false;
state.dtNext = state.dt;
while true
    inc = state.dt_increment;
    Y = Y0; it = 0; tol1 = Inf; tol2 = Inf; converged = false;
    while it < state.maxiter
        [Y, tol1, tol2, state] = OneIter(Y, Y0, state);
        it = it + 1;
        if ~isfinite(tol1), break; end
        if tol1 < toler, converged = true; break; end
        if state.cancelFlag, break; end
    end
    if converged || state.cancelFlag, break; end
    if state.dt <= dtMin, break; end              % give up; caller decides
    state.dt = max(dtMin, state.dt / inc);        % reject step, retry smaller
end

% Next-step dt from iteration count.
state.dtNext = state.dt;
if converged
    if it < 10
        state.dtNext = min(state.dt_max, state.dt * inc);
    elseif it > 20
        state.dtNext = state.dt / inc;
    end
end
state.dtNext = min(state.dtNext, state.dt_max); % enforce max dt on the next step
state.stepConverged = converged;
if ~converged
    Y = Y0;
    return;
end

state = refreshHydroFluxCache(state, Y);
state.Y = Y;
state.T = computeTemperatureProfile(Y, state);
end
