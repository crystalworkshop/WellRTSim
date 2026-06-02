function um = deriveMixtureVelocity(Y, state)
% Mass-weighted mixture velocity at each face i+1/2 from the 4-field state
% Y = [P; H; FVl; uv]:  u_mix = (rho_l*FVl + Sv*rho_v*uv) / rho_mix,
% with interface-averaged densities (matching RHS_v2 / convertToDriftFluxState).

n = size(Y, 2);
% [alpha_g, alpha_l, rho_g, rho_l, h_g, h_l, mu_l, T]
[Sv, ~, rv, rl, ~, ~, ~, ~] = calculatePhaseProperties(Y(1, :), Y(2, :), state);
rm = (1 - Sv).*rl + Sv.*rv;
FVl = Y(3, :); uv = Y(4, :);
um = zeros(1, n);
for i = 1:n
    iR = min(n, i+1);
    rl_f = 0.5*(rl(i) + rl(iR));
    Svrv_f = 0.5*(Sv(i)*rv(i) + Sv(iR)*rv(iR));
    rm_f = 0.5*(rm(i) + rm(iR));
    um(i) = (rl_f*FVl(i) + Svrv_f*uv(i)) / max(rm_f, 1e-9);
end
end
