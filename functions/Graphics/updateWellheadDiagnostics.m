function state = updateWellheadDiagnostics(state, updateFlowRate)
% Compute the latest wellhead sample from the drift-flux primaries at the
% wellhead face (state.Y(:,end) = [P; H; FVl; uv]) and store it in
% state.wellhead. The full time series lives only in the HDF5 file.

if nargin < 2 || isempty(updateFlowRate)
    updateFlowRate = true;
end

% [alpha_g, alpha_l, rho_g, rho_l, h_g, h_l, mu_l, T] at the top cell centre
[alpha_g, alpha_l, rho_g, rho_l, ~, ~, ~, Ttop] = ...
    calculatePhaseProperties(state.Y(1, end), state.Y(2, end), state);

FVl = state.Y(3, end);   % liquid volume flux at wellhead face
uv  = state.Y(4, end);   % vapour velocity at wellhead face
ul  = FVl / max(alpha_l, 1e-9);
ug  = uv;

areaTop = pi * state.Dp(end)^2 / 4;
rho_mix = alpha_g*rho_g + alpha_l*rho_l;
den = max(rho_mix, eps);

mdot_g = alpha_g * rho_g * ug * areaTop;
mdot_l = rho_l * FVl * areaTop;          % alpha_l*rho_l*ul = rho_l*FVl
denf = max(abs(mdot_g) + abs(mdot_l), eps);

state.wellhead.time_days = state.tt / 86400;
state.wellhead.WHP       = state.Y(1, end);
if updateFlowRate
    state.wellhead.flow_rate = mdot_g + mdot_l;
end
state.wellhead.rho_mix    = rho_mix;
state.wellhead.rho_gas    = rho_g;
state.wellhead.rho_liq    = rho_l;
state.wellhead.T          = Ttop;
state.wellhead.quality    = (alpha_g * rho_g) / den;
state.wellhead.steam_frac = abs(mdot_g) / denf;
state.wellhead.u_gas      = ug;
state.wellhead.u_liq      = ul;
state.wellhead.u_mix      = (rho_l*FVl + alpha_g*rho_g*uv) / den;
end
