function geom = chemistryPhaseGeometryV2(state)
% chemistryPhaseGeometryV2 Compute phase geometry/velocities for chemistry transport.

n = state.n;
P = state.Y(1, 1:n);
H = state.Y(2, 1:n);
U = state.Y(3, 1:n);

[alpha_g, alpha_l, rho_g, rho_l, ~, ~, ~, T] = calculatePhaseProperties(P, H, state);

Dp = state.Dp(:).';
Acell = pi * Dp.^2 / 4;
Vcell = Acell * state.dx;

Vl = alpha_l(:).' .* Vcell;
Vg = alpha_g(:).' .* Vcell;
VgL = Vg * 1000;
mLiquid = rho_l(:).' .* Vl;
mGas = rho_g(:).' .* Vg;

thetaNodes = state.gravityThetaNode(:).';
[u_g, u_l] = calculatePhaseVelocities(U, alpha_g, rho_g, rho_l, Dp, T, thetaNodes, state);

geom = struct();
geom.n = n;
geom.P_Pa = P(:).';
geom.P_bar = max(P(:).' / 1e5, 1e-8);
geom.T_K = T(:).';
geom.T_C = T(:).' - 273.15;
geom.alpha_g = alpha_g(:).';
geom.alpha_l = alpha_l(:).';
geom.rho_g = rho_g(:).';
geom.rho_l = rho_l(:).';
geom.u_g = u_g;
geom.u_l = u_l;
geom.Dp = Dp;
geom.Acell = Acell;
geom.Vcell = Vcell;
geom.Vl = Vl;
geom.Vg = Vg;
geom.Vg_L = VgL;
geom.mLiquid = mLiquid;
geom.mGas = mGas;
end
