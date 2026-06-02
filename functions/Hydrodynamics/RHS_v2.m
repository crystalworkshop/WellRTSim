function [S, state] = RHS_v2(Z3, Z2, Z1, Y0, state, i)
% 4-equation drift-flux residual (Tonkin et al. 2023 formulation).
%
% Primaries per block i:  Y(:,i) = [P_i; H_i; FVl_{i+1/2}; uv_{i+1/2}]
%   P,H : pressure [Pa], mixture enthalpy [J/kg] at the centre of cell i
%   FVl : liquid volume flux Sl*ul [m/s] at the face i+1/2 (right face of cell i)
%   uv  : vapour velocity [m/s] at the face i+1/2
%
% Residuals: S = [ mass(cell i); energy(cell i); momentum(i+1/2); slip(i+1/2) ].
% The slip equation closes the system implicitly so both phases' flow
% directions are known within each Newton step (no algebraic post-hoc slip).

S = zeros(4, 1);
avg = @(a,b) 0.5*(a+b);
dx = state.dx; dt = state.dt; g = state.g; n = state.n;

epsu = 1e-3;
if isfield(state, 'epsu_upwind') && ~isempty(state.epsu_upwind)
    epsu = state.epsu_upwind;
end

isBottom = (i == 1);
isTop = (i == n);
if isTop
    Z3 = buildTopGhostCell(state, Z2);
end

% ---- primaries ----
P2 = Z2(1); H2 = Z2(2); FVlR = Z2(3); uvR = Z2(4);   % cell i, face i+1/2
P1 = Z1(1); H1 = Z1(2); FVlL = Z1(3); uvL = Z1(4);   % cell i-1, face i-1/2
P3 = Z3(1); H3 = Z3(2); FVlRR = Z3(3); uvRR = Z3(4); % cell i+1, face i+3/2

% old-time values (cell i and right neighbour)
P2o = Y0(1,i); H2o = Y0(2,i); FVlRo = Y0(3,i); uvRo = Y0(4,i);
FVlLo = Y0(3,max(1,i-1)); uvLo = Y0(4,max(1,i-1));

% ---- geometry ----
iL = max(1, i-1); iR = min(n, i+1);
D2 = state.Dp(i); A2 = pi*D2^2/4;
A_im = pi*avg(state.Dp(iL), D2)^2/4;   % area at face i-1/2
A_ip = pi*avg(D2, state.Dp(iR))^2/4;   % area at face i+1/2
D_ip = avg(D2, state.Dp(iR));
A2o = pi*state.Dpprev(i)^2/4;
A_ip_o = pi*avg(state.Dpprev(i), state.Dpprev(iR))^2/4;

cth_im = state.gravityCthFace(i);
cth_ip = state.gravityCthFace(i+1);
th_ip = state.gravityThetaNode(i);

% ---- cell-centre thermodynamic state ----
% calculatePhaseProperties returns [alpha_g, alpha_l, rho_g, rho_l, h_g, h_l, mu_l, T]
[Sv2,Sl2,rv2,rl2,hv2,hl2,mu2,T2] = calculatePhaseProperties(P2,H2,state); rm2 = Sl2*rl2 + Sv2*rv2;
[Sv3,Sl3,rv3,rl3,hv3,hl3,mu3,T3] = calculatePhaseProperties(P3,H3,state); rm3 = Sl3*rl3 + Sv3*rv3;
[Sv1,Sl1,rv1,rl1,hv1,hl1,~,~] = calculatePhaseProperties(P1,H1,state);

% old cell-centre state (cell i and i+1) for transient terms
[Sv2o,Sl2o,rv2o,rl2o,hv2o,hl2o,~,~] = calculatePhaseProperties(P2o,H2o,state); rm2o = Sl2o*rl2o + Sv2o*rv2o;
P3o = Y0(1,iR); H3o = Y0(2,iR);
[Sv3o,Sl3o,rv3o,rl3o,~,~,~,~] = calculatePhaseProperties(P3o,H3o,state); rm3o = Sl3o*rl3o + Sv3o*rv3o;

% ---- feed ----
[Q_mass,Q_energy,state] = computeFeedSource(i, P2, A2, dx, state);

% =================== face i+1/2 fluxes (between cells i and i+1) ===================
Sl_ip = max(upw(FVlR, Sl2, Sl3, epsu), 1e-9);
rl_ip = upw(FVlR, rl2, rl3, epsu);
hl_ip = upw(FVlR, hl2, hl3, epsu);
Svrv_ip = upw(uvR, Sv2*rv2, Sv3*rv3, epsu);
hv_ip = upw(uvR, hv2, hv3, epsu);
ul_ip = FVlR / Sl_ip;
% Wellhead (M+1/2) production BC: the upwinded densities, saturations and
% enthalpies are ALWAYS taken from the block below (block M), regardless of
% the computed flow direction (Tonkin 2023, Section 4.5).
if isTop
    Sl_ip = max(Sl2, 1e-9);
    rl_ip = rl2; hl_ip = hl2;
    Svrv_ip = Sv2*rv2; hv_ip = hv2;
    ul_ip = FVlR / Sl_ip;
end
Fmv_ip = A_ip * Svrv_ip * uvR;     % vapour mass flux through face i+1/2
Fml_ip = A_ip * rl_ip * FVlR;      % liquid mass flux  (Sl*rl*ul = rl*FVl)
Fm_ip = Fmv_ip + Fml_ip;

% Wellhead no-flow switch (Tonkin 2023, Sections 4.5.3-4.5.5): if the wellhead
% total mass flux from the previous step is non-positive (production has
% stalled / would reverse), shut the well -- replace momentum with uv=0 and
% slip with FVl=0. Based on the old-time flux so it does not chatter within a
% Newton solve.
wellShut = false;
if isTop
    Fm_wh_old = A_ip*(Sv2o*rv2o*uvRo + rl2o*FVlRo);
    wellShut = (Fm_wh_old < 0);   % genuine backflow (not startup from rest)
end

% =================== face i-1/2 fluxes (between cells i-1 and i) ===================
Sl_im = max(upw(FVlL, Sl1, Sl2, epsu), 1e-9);
rl_im = upw(FVlL, rl1, rl2, epsu);
hl_im = upw(FVlL, hl1, hl2, epsu);
Svrv_im = upw(uvL, Sv1*rv1, Sv2*rv2, epsu);
hv_im = upw(uvL, hv1, hv2, epsu);
ul_im = FVlL / Sl_im;
Fmv_im = A_im * Svrv_im * uvL;
Fml_im = A_im * rl_im * FVlL;
if isBottom
    Fmv_im = 0; Fml_im = 0; ul_im = 0;
end
Fm_im = Fmv_im + Fml_im;

% record phase fluxes for diagnostics / consumers
state.Q_v_face(i) = Fmv_im; state.Q_l_face(i) = Fml_im;
state.Q_v_face(i+1) = Fmv_ip; state.Q_l_face(i+1) = Fml_ip;
state.Q_v(i) = 0.5*(Fmv_im + Fmv_ip); state.Q_l(i) = 0.5*(Fml_im + Fml_ip);
% per-phase upwind weights for the chemistry transport (weight on the
% left/lower cell), consistent with the upw() blend used above: each phase
% follows the sign of its own primary flux (FVl for liquid, uv for vapour),
% so near the boiling front the phases can upwind from opposite directions.
state.w_l_face(i)   = 0.5*(1 + tanh(FVlL/epsu));
state.w_v_face(i)   = 0.5*(1 + tanh(uvL/epsu));
state.w_l_face(i+1) = 0.5*(1 + tanh(FVlR/epsu));
state.w_v_face(i+1) = 0.5*(1 + tanh(uvR/epsu));
if isBottom, state.Q_v(i) = Fmv_ip; state.Q_l(i) = Fml_ip; end
if isTop,    state.Q_v(i) = Fmv_im; state.Q_l(i) = Fml_im; end

% =================== (1) mass residual, cell i ===================
S(1) = A2*rm2 - A2o*rm2o + dt*((Fm_ip - Fm_im)/dx - A2*Q_mass);

% =================== (2) energy residual, cell i ===================
% block-centred kinetic energy = average of squared face velocities (Eq. 61)
EKv2 = 0.25*(uvR^2 + uvL^2);
EKl2 = 0.25*(ul_ip^2 + ul_im^2);
e2  = Sv2*rv2*(hv2 + EKv2) + Sl2*rl2*(hl2 + EKl2) - P2;
ulRo = FVlRo / max(Sl2o, 1e-9); ulLo = FVlLo / max(Sl2o, 1e-9);
EKv2o = 0.25*(uvRo^2 + uvLo^2);
EKl2o = 0.25*(ulRo^2 + ulLo^2);
e2o = Sv2o*rv2o*(hv2o + EKv2o) + Sl2o*rl2o*(hl2o + EKl2o) - P2o;
Etr = A2*e2 - A2o*e2o;

jh_ip = Fmv_ip*(hv_ip + 0.5*uvR^2) + Fml_ip*(hl_ip + 0.5*ul_ip^2);
jh_im = Fmv_im*(hv_im + 0.5*uvL^2) + Fml_im*(hl_im + 0.5*ul_im^2);
grav_work = g*avg(cth_im*Fm_im, cth_ip*Fm_ip);
q_heat = heatLoss(i, A2, D2, T2, state); state.q_heat(i) = q_heat;
state.Energy(i) = e2; state.Energy_flux(i) = (jh_ip - jh_im)/dx;
S(2) = Etr + dt*((jh_ip - jh_im)/dx + grav_work + q_heat - A2*Q_energy);

% =================== (3) momentum residual, face i+1/2 ===================
% mixture momentum per area at the face = Sl*rl*ul + Sv*rv*uv = rl*FVl + Svrv*uv
mom_ip = rl_ip*FVlR + Svrv_ip*uvR;
Svrv_ip_o = upw(uvRo, Sv2o*rv2o, Sv3o*rv3o, epsu);
rl_ip_o   = upw(FVlRo, rl2o, rl3o, epsu);
mom_ip_o  = rl_ip_o*FVlRo + Svrv_ip_o*uvRo;

% Convective momentum flux difference, one-sided per flow direction (Eqs. 77-79):
% up-flow uses the backward difference (i+1/2)-(i-1/2); down-flow the forward
% difference (i+3/2)-(i+1/2). This keeps the momentum advection upwind (no
% central-difference oscillations) rather than the cell-centred blend.
iRR = min(n, i+2);
A_ipp = pi*avg(state.Dp(iR), state.Dp(iRR))^2/4;
ul_ipp = FVlRR / max(Sl3, 1e-9);   % liquid velocity at face i+3/2 (Sl from cell i+1)
Slrl_im = upw(FVlL, Sl1*rl1, Sl2*rl2, epsu);
Slrl_ip = upw(FVlR, Sl2*rl2, Sl3*rl3, epsu);
% momentum flux M = A*(S*rho)*u^2 at velocity nodes i-1/2, i+1/2, i+3/2
Mim_g = A_im*Svrv_im*uvL^2;    Mip_g = A_ip*Svrv_ip*uvR^2;    Mipp_g = A_ipp*Sv3*rv3*uvRR^2;
Mim_l = A_im*Slrl_im*ul_im^2;  Mip_l = A_ip*Slrl_ip*ul_ip^2;  Mipp_l = A_ipp*Sl3*rl3*ul_ipp^2;
if isBottom
    Mim_g = 0; Mim_l = 0;   % closed bottom: no momentum flux through face 1/2
end
wmg = 0.5*(1 + tanh(uvR/epsu));
wml = 0.5*(1 + tanh(ul_ip/epsu));
dFmom_g = wmg*(Mip_g - Mim_g) + (1 - wmg)*(Mipp_g - Mip_g);
dFmom_l = wml*(Mip_l - Mim_l) + (1 - wml)*(Mipp_l - Mip_l);
dFmom = dFmom_g + dFmom_l;

rho_e = avg(rm2, rm3);
umix_ip = mom_ip / max(rho_e, 1e-9);
mu_e = avg(mu2, mu3);
Re_e = max(1e-1, rho_e*abs(umix_ip)*D_ip/max(mu_e, 1e-12));
rough_e = avg(state.eps(i), state.eps(iR));
fD = frictionFactor(D_ip, Re_e, rough_e);
S_fric = A_ip * fD * rho_e * umix_ip*abs(umix_ip) / (2*max(D_ip, 1e-12));
S_grav = rho_e * g * A_ip * cth_ip;
dp = A_ip*(P3 - P2);
if isTop && wellShut
    S(3) = uvR;   % shut-in wellhead: no vapour flow (Tonkin 2023, Eq 92)
else
    S(3) = A_ip*mom_ip - A_ip_o*mom_ip_o + dt*((dp + dFmom)/dx + S_fric + S_grav);
end

% =================== (4) slip residual, face i+1/2 ===================
% uv = FV + ud  with  FV = FVl + Sv*uv  (C0 = 1) -> uv*(1-Sv) - FVl - ud = 0.
% Drift velocity via the Pan/Shi interface blend (Eqs. 81-88): udS uses the
% saturation below the face, udN the saturation above, blended by Uslip(Save)
% over Save in [0.3,0.4]. ud thus depends on Sv_i and Sv_{i+1} individually
% (damping odd-even Sv oscillations) while staying smooth in velocity.
rl_f = avg(rl2, rl3); rv_f = avg(rv2, rv3); T_f = avg(T2, T3);
Save = avg(Sv2, Sv3);                    % average Sv: Uslip blend only (Eq 87)
Sv_up = upw(uvR, Sv2, Sv3, epsu);        % upstream-weighted interface Sv (Eq 80)
udS = driftVelocity(Sv2, rl_f, rv_f, T_f, D_ip, th_ip, state);
udN = driftVelocity(Sv3, rl_f, rv_f, T_f, D_ip, th_ip, state);
a1s = 0.3; a2s = 0.4;
if Save <= a1s
    Uslip = 0;
elseif Save < a2s
    Uslip = 0.5*(1 - cos(pi*(Save - a1s)/(a2s - a1s)));
else
    Uslip = 1;
end
ud = (1 - Uslip)*udS + Uslip*udN;
% Slip residual (Eq 80): the (Sv*uv) term uses the UPWINDED interface Sv, not
% the central average -- a central Sv here decouples odd/even cells and
% produces a checkerboard oscillation in the velocity field.
S_drift = uvR*(1 - Sv_up) - FVlR - ud;

% Disappearing-liquid limit (Save -> 1): the drift form's uv coefficient
% (1-Save) -> 0, so uv loses its constraint and the Jacobian goes singular.
% Blend to a no-slip closure ul = uv, written in volume-flux form
% FVl = Sl*uv, which cleanly pins the vanishing liquid flux (coefficient 1 on
% FVl) instead of leaving uv undetermined. Active only for Save > 0.98.
% a1n = 0.98; a2n = 0.999;
% if Save <= a1n
%     wNoSlip = 0;
% elseif Save < a2n
%     wNoSlip = 0.5*(1 - cos(pi*(Save - a1n)/(a2n - a1n)));
% else
%     wNoSlip = 1;
% end
wNoSlip = 0;
S_noslip = FVlR - Sl_ip*uvR;
S(4) = (1 - wNoSlip)*S_drift + wNoSlip*S_noslip;

% Wellhead slip BC (Tonkin 2023, Eqs 80 & 90): a simple drift velocity from
% block M with the saturation from below, NOT the interior two-sided blend.
if isTop
    ud_M = driftVelocity(Sv2, rl2, rv2, T2, D_ip, th_ip, state);
    S(4) = uvR*(1 - Sv2) - FVlR - ud_M;
end
% Shut-in override (Eq 93): no liquid flow through a stalled wellhead.
if isTop && wellShut
    S(4) = FVlR;
end

end

% ---------------------------------------------------------------------------
function v = upw(dir, vL, vR, epsu)
% Smoothed upwind: dir>0 takes the left (upstream) value.
w = 0.5*(1 + tanh(dir/epsu));
v = w*vL + (1-w)*vR;
end

function F = FLUX(VL, VR, U, epsu)
w = 0.5*(1 + tanh(U/epsu));
F = (w*VL + (1-w)*VR)*U;
end

function q_heat = heatLoss(i, A2, D2, T2, state)
if state.tt <= state.t_adjust + state.t_transit/2
    q_heat = 0;
    return;
end
kappa = state.k_r / max(state.C_r*state.rho_r, 1e-12);
tD = kappa*state.tt / max((D2/2)^2, 1e-12);
ftd = 0.982*log(1 + 1.81*tD);
iRock = min(max(i, 1), numel(state.T_rock));
H_q = state.H_q_node(i);
q_heat = A2*4/max(D2, 1e-12) * H_q*state.k_r*(T2 - state.T_rock(iRock)) / ...
    max(state.k_r + H_q*D2*ftd/2, 1e-12);
end

function fD = frictionFactor(D, Re, roughness)
Re = max(Re, 1e-12);
roughness = max(roughness, 0);
term = (2e4 * roughness / max(D, 1e-12) + 1e6 / Re)^(1/3);
fD = 0.0055 * (1 + term);
end

function Zg = buildTopGhostCell(state, Z2)
% Pressure outlet boundary condition (iBC_top = 3): prescribe only the
% wellhead pressure. Enthalpy, liquid volume flux and vapour velocity are
% zero-gradient (taken from the top cell), so for outflow the wellhead face
% sees the interior steam-rich state instead of a spurious liquid ghost,
% which previously collapsed the outlet velocity.
Zg = [state.P_top; Z2(2); Z2(3); Z2(4)];
end

function [Q_mass,Q_energy,state] = computeFeedSource(i, pCell, ACell, dx, state)
Q_mass = 0;
Q_energy = 0;
state.Q_mass(i) = 0;
if state.feedzone_cells(i) == 1
    if state.feed == 3
        s_m = state.feedzone_Qm(i);
    else
        s_m = state.feedzone_PI(i) * (state.feedzone_P_res(i) - pCell);
    end
    Q_mass = s_m / (ACell*dx);
    Q_energy = state.feedzone_H_res(i) * Q_mass;
    state.Q_mass(i) = Q_mass;
end
end
