function Y4 = convertToDriftFluxState(state, Y3)
% Convert a 3-field steady profile Y3 = [P; H; u_mix] (u_mix at faces i+1/2)
% into the 4-field drift-flux layout Y4 = [P; H; FVl; uv].
%
% Uses the slip relation (C0 = 1) so the slip residual is satisfied exactly:
%   uv  = u_mix + (rho_l/rho_mix) * ud
%   FVl = uv*(1 - Sv) - ud
% with interface-averaged Sv, rho and the Shi drift velocity ud, matching the
% face quantities used in RHS_v2.

n = state.n;
P = Y3(1, :); H = Y3(2, :); um = Y3(3, :);
Y4 = zeros(4, n);
Y4(1, :) = P;
Y4(2, :) = H;

% [alpha_g, alpha_l, rho_g, rho_l, h_g, h_l, mu_l, T]
[Sv, ~, rv, rl, ~, ~, ~, T] = calculatePhaseProperties(P, H, state);
rmix = (1 - Sv).*rl + Sv.*rv;

for i = 1:n
    iR = min(n, i+1);
    Sv_f = 0.5*(Sv(i) + Sv(iR));
    rl_f = 0.5*(rl(i) + rl(iR));
    rv_f = 0.5*(rv(i) + rv(iR));
    rm_f = 0.5*(rmix(i) + rmix(iR));
    T_f  = 0.5*(T(i) + T(iR));
    D_ip = 0.5*(state.Dp(i) + state.Dp(iR));
    th_ip = state.gravityThetaNode(i);

    ud = driftVelocity(Sv_f, rl_f, rv_f, T_f, D_ip, th_ip, state);
    uv = um(i) + (rl_f/max(rm_f,1e-9))*ud;
    FVl = uv*(1 - Sv_f) - ud;
    Y4(3, i) = FVl;
    Y4(4, i) = uv;
end
end
