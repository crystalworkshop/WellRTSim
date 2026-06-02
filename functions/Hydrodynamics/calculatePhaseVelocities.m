function [u_g, u_l] = calculatePhaseVelocities(u_mix, alpha_g, rho_v, rho_l, d, T, theta0, state)
u_g = u_mix;
u_l = u_mix;
% return; % Disable for now until we can validate against experimental data
active = state.tt >= state.t_adjust + state.t_transit & alpha_g > 1e-5;
if ~any(active(:))
    return;
end

eps_h = 1e-6;
Sv = min(1 - eps_h, max(eps_h, alpha_g));
% IAPWS-2014 surface tension of ordinary water: sigma = B*tau^mu*(1 - b*tau),
% B = 0.2358 N/m, mu = 1.256, b = 0.625, tau = 1 - T/Tc, Tc = 647.096 K.
sigma = zeros(size(T));
Tc = 647.096;
subcritical = T < Tc;
tau = 1 - T(subcritical) / Tc;
sigma(subcritical) = 0.2358 * tau.^1.256 .* (1 - 0.625 * tau);
sigma = max(sigma, 1e-8);

alpha1 = 0.06;
alpha2 = 0.21;
Cku = 142;
Cw = 0.008;
g = state.g;

deltaRho = max(rho_l - rho_v, 0);
uc = (g * sigma .* deltaRho ./ max(rho_l.^2, 1e-12)).^(1/4);
Nb = d.^2 * g .* deltaRho ./ sigma;
Ku = sqrt((Cku ./ max(sqrt(Nb), 1e-12)) .* (sqrt(1 + Nb / (Cku^2 * Cw)) - 1));

K = Ku;
lowHoldup = Sv <= alpha1;
midHoldup = Sv > alpha1 & Sv < alpha2;
K(lowHoldup) = 1.53;
K(midHoldup) = 1.53 + 0.5 * (Ku(midHoldup) - 1.53) .* ...
    (1 - cos(pi * (Sv(midHoldup) - alpha1) / (alpha2 - alpha1)));

m_theta = 1.85 * (cos(theta0).^0.21) .* (1 + sin(theta0)).^0.95;
den = Sv .* sqrt(max(rho_v ./ max(rho_l, 1e-12), 0)) + (1 - Sv);
ud = m_theta .* uc .* ((1 - Sv) .* K ./ max(den, eps_h));

rho_mix = (1 - Sv) .* rho_l + Sv .* rho_v;
Fv = u_mix + Sv .* (rho_l - rho_v) ./ max(rho_mix, eps_h) .* ud;

fullGas = active & Sv >= 1 - eps_h;
twoPhase = active & Sv > eps_h & ~fullGas;

u_g(fullGas) = u_mix(fullGas) + ud(fullGas);
u_g(twoPhase) = Fv(twoPhase) + ud(twoPhase);
u_l(twoPhase) = Fv(twoPhase) - (Sv(twoPhase) ./ max(1 - Sv(twoPhase), eps_h)) .* ud(twoPhase);
end
