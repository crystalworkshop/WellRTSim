function [S,state] = RHS_v2(Z3, Z2, Z1, Y0, state, i)
S = zeros(3,1);
avg = @(a,b) 0.5*(a+b);

dx = state.dx;
dt = state.dt;
g = state.g;
n = state.n;
epsu = 1e-6;

isBottom = (i == 1);
isTop = (i == n);

state.Q_mass(i) = 0;

if isBottom
    Z1 = [Z2(1); Z2(2); 0];
end
if isTop
    Z3 = buildTopGhostCell(state, Y0(:,n));
end

p1 = Z1(1); h1 = Z1(2); u1 = Z1(3);
p2 = Z2(1); h2 = Z2(2); u2 = Z2(3);
p3 = Z3(1); h3 = Z3(2); u3 = Z3(3);

p20 = Y0(1,i);
h20 = Y0(2,i);
u20 = Y0(3,i);
iL = max(1, i-1);
p10 = Y0(1,iL);
h10 = Y0(2,iL);

D1 = state.Dp(iL);
D2 = state.Dp(i);
A2 = pi*D2^2/4;
iR = min(n, i+1);
D3 = state.Dp(iR);

D20 = state.Dpprev(i);
A20 = pi*D20^2/4;
D120 = avg(state.Dpprev(iL), state.Dpprev(i));
A120 = pi*D120^2/4;
D12 = avg(D1, D2);
A12 = pi*D12^2/4;
D23 = avg(D2, D3);
A23 = pi*D23^2/4;

theta1 = state.gravityThetaNode(iL);
theta2 = state.gravityThetaNode(i);
theta3 = state.gravityThetaNode(iR);
cth12 = state.gravityCthFace(i);
cth23 = state.gravityCthFace(i + 1);

[ag1,al1,rg1,rl1,hg1,hl1,mu1,T1] = calculatePhaseProperties(p1,h1,state);
[ug1,ul1] = calculatePhaseVelocities(u1,ag1,rg1,rl1,D12,T1,theta1,state);
rho1 = ag1*rg1 + al1*rl1;

[ag2,al2,rg2,rl2,hg2,hl2,mu2,T2] = calculatePhaseProperties(p2,h2,state);
[ug2,ul2] = calculatePhaseVelocities(u2,ag2,rg2,rl2,D12,T2,theta2,state);
rho2 = ag2*rg2 + al2*rl2;

[ag3,al3,rg3,rl3,hg3,hl3,~,T3] = calculatePhaseProperties(p3,h3,state);
[ug3,ul3] = calculatePhaseVelocities(u3,ag3,rg3,rl3,D23,T3,theta3,state);

[ag20,al20,rg20,rl20,hg20,hl20,~,T20] = calculatePhaseProperties(p20,h20,state);
[ug20,ul20] = calculatePhaseVelocities(u20,ag20,rg20,rl20,D12,T20,theta2,state);
rho20 = ag20*rg20 + al20*rl20;

[ag10,al10,rg10,rl10] = calculatePhaseProperties(p10,h10,state);
rho10 = ag10*rg10 + al10*rl10;

[Q_mass,Q_energy,state] = computeFeedSource(i, p2, A2, dx, state);

w12g = 0.5 * (1 + tanh(ug1 / epsu));
w12l = 0.5 * (1 + tanh(ul1 / epsu));
w23g = 0.5 * (1 + tanh(ug2 / epsu));
w23l = 0.5 * (1 + tanh(ul2 / epsu));

mdot12_g = FLUX(A12*ag1*rg1, A12*ag2*rg2, ug1, epsu);
mdot12_l = FLUX(A12*al1*rl1, A12*al2*rl2, ul1, epsu);
if isBottom
    mdot12_g = 0;
    mdot12_l = 0;
end
mdot12 = mdot12_g + mdot12_l;

mdot23_g = FLUX(A23*ag2*rg2, A23*ag3*rg3, ug2, epsu);
mdot23_l = FLUX(A23*al2*rl2, A23*al3*rl3, ul2, epsu);
mdot23 = mdot23_g + mdot23_l;

state.Q_v_face(i) = mdot12_g;
state.Q_l_face(i) = mdot12_l;
state.w_v_face(i) = w12g;
state.w_l_face(i) = w12l;
state.Q_v_face(i + 1) = mdot23_g;
state.Q_l_face(i + 1) = mdot23_l;
state.w_v_face(i + 1) = w23g;
state.w_l_face(i + 1) = w23l;

state.Q_v(i) = 0.5*(mdot12_g + mdot23_g);
state.Q_l(i) = 0.5*(mdot12_l + mdot23_l);
if isBottom, state.Q_v(i) = mdot23_g; state.Q_l(i) = mdot23_l; end
if isTop, state.Q_v(i) = mdot12_g; state.Q_l(i) = mdot12_l; end

mdiv = (mdot23 - mdot12) / dx;
S(1) = A2*rho2 - A20*rho20 + dt*(mdiv - A2*Q_mass);

ug_12 = avg(ug1, ug2);
ul_12 = avg(ul1, ul2);
ug_23 = avg(ug3, ug2);
ul_23 = avg(ul3, ul2);

Fmom12_g = FLUX(A12*ag1*rg1*ug1, A12*ag2*rg2*ug2, ug_12, epsu);
Fmom12_l = FLUX(A12*al1*rl1*ul1, A12*al2*rl2*ul2, ul_12, epsu);
if isBottom
    Fmom12_g = 0;
    Fmom12_l = 0;
end
Fmom12 = Fmom12_g + Fmom12_l;

Fmom23_g = FLUX(A23*ag2*rg2*ug2, A23*ag3*rg3*ug3, ug_23, epsu);
Fmom23_l = FLUX(A23*al2*rl2*ul2, A23*al3*rl3*ul3, ul_23, epsu);
Fmom23 = Fmom23_g + Fmom23_l;

rho_e = avg(rho1, rho2);
mu_e = avg(mu1, mu2);
u_e = u2;
Re_e = max(1e-1, rho_e*abs(u_e)*D12/max(mu_e, 1e-12));
rough_e = state.eps(i);
fD = frictionFactor(D12, Re_e, rough_e);
S_fric = A12 * fD * rho_e * u_e*abs(u_e) / (2*max(D12, 1e-12));
S_grav = rho_e * g * A12 * cth12;

dp = A12*(p3 - p2);
rho_e0 = avg(rho10, rho20);
S(2) = A12*rho_e*u2 - A120*rho_e0*u20 + ...
    dt*((dp + Fmom23 - Fmom12)/dx + S_fric + S_grav);

e2 = ag2*rg2*(hg2 + 0.5*ug2^2) + al2*rl2*(hl2 + 0.5*ul2^2) - p2;
e20 = ag20*rg20*(hg20 + 0.5*ug20^2) + al20*rl20*(hl20 + 0.5*ul20^2) - p20;
Etr = A2*e2 - A20*e20;

jh12 = FLUX(hg1 + 0.5*ug2^2, hg2 + 0.5*ug2^2, mdot12_g, epsu) + ...
    FLUX(hl1 + 0.5*ul2^2, hl2 + 0.5*ul2^2, mdot12_l, epsu);
if isBottom
    jh12 = 0;
end
jh23 = FLUX(hg2 + 0.5*ug3^2, hg3 + 0.5*ug3^2, mdot23_g, epsu) + ...
    FLUX(hl2 + 0.5*ul3^2, hl3 + 0.5*ul3^2, mdot23_l, epsu);

grav_work = g*avg(cth12*mdot12, cth23*mdot23);

if state.tt <= state.t_adjust + state.t_transit/2
    state.q_heat(i) = 0;
else
    kappa = state.k_r / max(state.C_r*state.rho_r, 1e-12);
    tD = kappa*state.tt / max(D2^2, 1e-12);
    ftd = 0.982*log(1 + 1.81*tD);
    iRock = min(max(i, 1), numel(state.T_rock));
    H_q = state.H_q_node(i);
    state.q_heat(i) = A2*4/max(D2, 1e-12) * H_q*state.k_r*(T2 - state.T_rock(iRock))/  ...
        max(state.k_r + H_q*D2*ftd/2, 1e-12);
    
end

state.Energy(i) = e2;
state.Energy_flux(i) = (jh23 - jh12)/dx;
S(3) = Etr + dt*((jh23 - jh12)/dx + grav_work + state.q_heat(i) - A2*Q_energy);

if any(isnan(S))
    disp(num2str(Z1', '%.16g\t'));
    disp(num2str(Z2', '%.16g\t'));
    disp(num2str(Z3', '%.16g\t'));
end
end

function F = FLUX(VL, VR, U, epsu)
w = 0.5*(1 + tanh(U/epsu));
F = (w*VL + (1-w)*VR)*U;
end

function fD = frictionFactor(D, Re, roughness)
Re = max(Re, 1e-12);
roughness = max(roughness, 0);
term = (2e4 * roughness / max(D, 1e-12) + 1e6 / Re)^(1/3);
fD = 0.0055 * (1 + term);
end

function Zg = buildTopGhostCell(state, yTopPrev)
pGhost = state.P_top;
hGhost = yTopPrev(2);
uGhost = yTopPrev(3);

if state.H_top > 0
    if state.H_top < 1e5
        hGhost = 1e3*state.H_top;
    else
        hGhost = state.H_top;
    end
end

Zg = [pGhost; hGhost; uGhost];
end

function [Q_mass,Q_energy,state] = computeFeedSource(i, pCell, ACell, dx, state)
Q_mass = 0;
Q_energy = 0;
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
