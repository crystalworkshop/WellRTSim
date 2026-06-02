function state = updateTopPressureFromWHP(state)
% Update top pressure according to top boundary-condition mode, then apply
% optional wellhead flow control.
%
% iBC_top modes:
%   1 - pressure from WHP schedule (legacy behavior)
%   2 - flow-rate mode (pressure not overridden here)
%   3 - fixed outlet pressure equal to P_top from setup
%
% Wellhead flow control (whp_flow_control == 1): if the most recent wellhead
% mass flow is negative (fluid flowing back down into the well), the wellhead
% pressure is shed toward atmospheric to restore positive production. The
% relief is one-way (WHP only drops) and floored at state.P_atm.

bcMode = 1;
if isfield(state, 'iBC_top') && isnumeric(state.iBC_top) && isfinite(state.iBC_top)
    bcMode = round(state.iBC_top);
end

% Flow-rate BC: pressure controller is handled elsewhere; nothing to do.
if bcMode == 2
    return;
end

% ---- base wellhead pressure for this mode ----
if bcMode == 3
    if ~isfield(state, 'P_top_fixed') || ~isnumeric(state.P_top_fixed) || ...
            ~isfinite(state.P_top_fixed)
        state.P_top_fixed = state.P_top;
    end
    pTop = state.P_top_fixed;
else
    pTop = whpScheduleValue(state);
end

% ---- optional wellhead flow control ----
if isfield(state, 'whp_flow_control') && state.whp_flow_control == 1
    Patm = state.P_atm;
    flow = NaN;
    if isfield(state, 'wellhead') && isfield(state.wellhead, 'flow_rate')
        flow = state.wellhead.flow_rate;
    end
    if isfinite(flow) && flow < 0
        % Latch the controlled ceiling to the current target on first backflow,
        % then shed a fraction of the gap to atmospheric each step.
        if ~isfinite(state.WHP_ctrl)
            state.WHP_ctrl = pTop;
        end
        state.WHP_ctrl = max(Patm, ...
            state.WHP_ctrl - state.whp_ctrl_relax * (state.WHP_ctrl - Patm));
    end
    if isfinite(state.WHP_ctrl)
        pTop = min(pTop, state.WHP_ctrl);
    end
    pTop = max(Patm, pTop);
end

state.P_top = pTop;
end

function pTop = whpScheduleValue(state)
% Wellhead pressure from the WHP schedule with soft-start (legacy behavior).
sched = state.whpSchedule;
t_adjust = state.t_adjust;
t_transit = state.t_transit;
t_eval = state.tt + state.dt;

if t_eval <= t_adjust
    pTop = sched.p_init;
    return;
end

p_first = sched.pPa(1);
if t_eval <= t_adjust + t_transit
    denom = max(t_transit, eps);
    frac = min(1, max(0, (t_eval - t_adjust) / denom));
    pTop = sched.p_init + frac * (p_first - sched.p_init);
    return;
end

t_rel = t_eval - t_adjust - t_transit;
if t_rel <= sched.tSec(1)
    pTop = p_first;
elseif t_rel >= sched.tSec(end)
    pTop = sched.pPa(end);
else
    pTop = interp1(sched.tSec, sched.pPa, t_rel, 'linear');
end
if ~isfinite(pTop)
    pTop = state.P_top;
end
end
