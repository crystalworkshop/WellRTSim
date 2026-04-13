function state = updateWellheadDiagnostics(state, idx, phaseVelocity, mixVelocity, updateFlowRate)
% Update cached wellhead time-series diagnostics for the current state.

if nargin < 4 || isempty(mixVelocity)
    mixVelocity = phaseVelocity;
end
if nargin < 5
    updateFlowRate = true;
end

state.tsav(idx) = state.tt / 86400;
state.WHP(idx) = state.Y(1, end);

[alpha_g, alpha_l, rho_g, rho_l, ~, ~, ~, Ttop] = ...
    calculatePhaseProperties(state.Y(1, end), state.Y(2, end), state);
thetaTop = state.gravityThetaNode(end);
[ug, ul] = calculatePhaseVelocities(phaseVelocity, alpha_g, rho_g, rho_l, ...
    state.Dp(end), Ttop, thetaTop, state);

areaTop = pi * state.Dp(end)^2 / 4;
if updateFlowRate
    state.FlowRate(idx) = ((1 - alpha_g) * rho_l * ul + alpha_g * rho_g * ug) * areaTop;
end

den = alpha_g * rho_g + alpha_l * rho_l;
if den == 0
    den = eps;
end

state.RhoTop(idx) = alpha_g * rho_g + alpha_l * rho_l;
state.RhoGasTop(idx) = rho_g;
state.RhoLiqTop(idx) = rho_l;
state.TTop(idx) = Ttop;
state.Quality(idx) = (alpha_g * rho_g) / den;

mdot_g = alpha_g * rho_g * ug * areaTop;
mdot_l = alpha_l * rho_l * ul * areaTop;
denf = abs(mdot_g) + abs(mdot_l);
if denf == 0
    denf = eps;
end

state.SteamFrac(idx) = abs(mdot_g) / denf;
state.UGasTop(idx) = ug;
state.ULiqTop(idx) = ul;
state.UMixTop(idx) = mixVelocity;
end
