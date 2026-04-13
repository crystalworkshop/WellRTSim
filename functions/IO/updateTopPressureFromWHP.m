function state = updateTopPressureFromWHP(state)
% Update top pressure according to top boundary-condition mode.
%
% iBC_top modes:
%   1 - pressure from WHP schedule (legacy behavior)
%   2 - flow-rate mode (pressure not overridden here)
%   3 - fixed outlet pressure equal to P_top from setup

sched = state.whpSchedule;

t_adjust = state.t_adjust;
t_transit = state.t_transit;

bcMode = 1;
if isfield(state, 'iBC_top') && isnumeric(state.iBC_top) && isfinite(state.iBC_top)
    bcMode = round(state.iBC_top);
end

if bcMode == 3
    % Keep pressure fixed at the configured outlet value.
    if ~isfield(state, 'P_top_fixed') || ~isnumeric(state.P_top_fixed) || ...
            ~isfinite(state.P_top_fixed)
        state.P_top_fixed = state.P_top;
    end
    state.P_top = state.P_top_fixed;
    return;
end

if bcMode == 2
    % Flow-rate BC: pressure controller is handled elsewhere.
    return;
end

t_eval = state.tt + state.dt;

if t_eval <= t_adjust
    state.P_top = sched.p_init;
    return;
end

p_first = sched.pPa(1);
if t_eval <= t_adjust + t_transit
    denom = max(t_transit, eps);
    frac = (t_eval - t_adjust) / denom;
    frac = min(1, max(0, frac));
    state.P_top = sched.p_init + frac * (p_first - sched.p_init);
    return;
end

t_rel = t_eval - t_adjust - t_transit;
if t_rel <= sched.tSec(1)
    p_target = p_first;
elseif t_rel >= sched.tSec(end)
    p_target = sched.pPa(end);
else
    p_target = interp1(sched.tSec, sched.pPa, t_rel, 'linear');
end

if isfinite(p_target)
    state.P_top = p_target;
end
end
