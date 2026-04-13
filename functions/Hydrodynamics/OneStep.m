function [Y, tol1, tol2, it, state] = OneStep(state)
toler = state.epsQ;
it = 0;
Y = state.Y0;
Y0 = state.Y0;
tol1 = Inf;
tol2 = Inf;

while it <= state.maxiter
    [Y, tol1, tol2, state] = OneIter(Y, Y0, state);
    if tol1 < toler || isnan(tol1) || state.dt < 1e-5 || state.cancelFlag
        break;
    end
    it = it + 1;
end

if it < 8
    state.dt = min(state.dt_max, state.dt * state.dt_increment);
elseif it > state.maxiter + 10
    state.dt = state.dt / state.dt_increment;
end

state = refreshHydroFluxCache(state, Y);
state.Y = Y;
state.T = computeTemperatureProfile(Y, state);
end
