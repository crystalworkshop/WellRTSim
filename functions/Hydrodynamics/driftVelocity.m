function ud = driftVelocity(Sv, rl, rv, T, d, theta, state)
% Shi et al. (2005) drift velocity (C0 = 1), evaluated at an interface.
%   Sv      vapour saturation (void fraction) at the interface
%   rl,rv   liquid / vapour densities [kg/m^3]
%   T       temperature [K]
%   d       pipe diameter [m]
%   theta   deviation from vertical [rad]
% Returns the drift velocity ud [m/s].

eps_h = 1e-6;
Sv = min(1 - eps_h, max(eps_h, Sv));

Tc = 647.096;                       % IAPWS-2014 surface tension
if T < Tc
    tau = 1 - T/Tc;
    sigma = 0.2358 * tau^1.256 * (1 - 0.625*tau);
else
    sigma = 1e-8;
end
sigma = max(sigma, 1e-8);

g = state.g; Cku = 142; Cw = 0.008; a1 = 0.06; a2 = 0.21;
dR = max(rl - rv, 0);
uc = (g*sigma*dR / max(rl^2, 1e-12))^(1/4);          % characteristic velocity (Eq. 32)
Nb = d^2 * g * dR / sigma;                            % Bond number (Eq. 35)
Ku = sqrt((Cku / max(sqrt(Nb), 1e-12)) * (sqrt(1 + Nb/(Cku^2*Cw)) - 1)); % Kutateladze (Eq. 34)

if Sv <= a1
    K = 1.53;
elseif Sv < a2
    K = 1.53 + 0.5*(Ku - 1.53)*(1 - cos(pi*(Sv - a1)/(a2 - a1)));
else
    K = Ku;
end

m = 1.85 * (cos(theta)^0.21) * (1 + sin(theta))^0.95;  % inclination factor (Eq. 36)
den = Sv*sqrt(max(rv/max(rl,1e-12), 0)) + (1 - Sv);
ud = m * uc * ((1 - Sv) * K / max(den, eps_h));        % Eq. 42 (C0 = 1)
end
