function state = initializeSteadyState(state)
% Steady-state initializer (clean version).
% IC_switch:
%   1 - top-down integration
%   2 - bottom-up shooting integration
%   3 - prescribed profiles: P from InitPressure.csv, T from
%       InitTemperature.csv (-> H via IAPWS), velocity from Q_init

model = buildSteadyModel(state);

switch state.IC_switch
    case 1
        [P, H, Mdot, ~] = solveTopDown(state, model);
    case 2
        [P, H, Mdot, pBotUsed] = solveBottomUp(state, model);
    case 3
        [P, H, Mdot, ~] = solvePrescribed(state, model);
    otherwise
        error('initializeSteadyState:InvalidICSwitch', ...
            'IC_switch must be 1 (top-down), 2 (bottom-up), or 3 (prescribed profiles).');
end

n = state.n;
U = zeros(1, n);
for i = 1:n
    U(i) = velocityFromMdot(Mdot(i), P(i), H(i), state.Dp(i), state, model.mCut);
end
Y = [P; H; U];

state.Y = Y;
state.Y0 = Y;
state.T = computeTemperatureProfile(Y, state);
state.H_bot = H(1);
state.P_top = P(end);
state.Q = Mdot(end);
state.P_res = P;
state.H_res = H;
state.Q_mass = zeros(1, n);
for i = 1:n
    s = state.Lp - state.x(i);
    qMassLine = sourceMassLine(s, P(i), model);
    Ai = pi * state.Dp(i)^2 / 4;
    state.Q_mass(i) = qMassLine / max(Ai, 1e-12);
end
state.Q_v = zeros(1, n);
state.Q_l = zeros(1, n);
state.mdot_ic = Mdot;

% Convert the 3-field steady profile [P;H;u_mix] to the drift-flux layout
% [P;H;FVl;uv] used by the transient solver (slip residual zero by construction).
Y = convertToDriftFluxState(state, Y);
state.k = 4;
state.Y = Y;
state.Y0 = Y;

% Populate phase-flow caches for the converged initial state so the
% initial profile plots show liquid/gas/mixture rates immediately.
state = refreshHydroFluxCache(state, Y);
state = updateWellheadDiagnostics(state, false);

Ttop = state.T(end) - 273.15;
pScale = state.pressureUnitScale;
pLabel = char(state.pressureUnitLabel);
fprintf('  Top:    P=%.3f %s, T=%.1f C, Q=%.3f kg/s\n', ...
    state.P_top/pScale, pLabel, Ttop, Mdot(end));
if state.IC_switch == 2 && pBotUsed > state.P_bot
    fprintf('  Pbot increased from %.3f %s to %.3f %s to reach top with P>=P_top.\n', ...
        state.P_bot/pScale, pLabel, pBotUsed/pScale, pLabel);
end
end

function model = buildSteadyModel(state)
model = struct();
model.state = state;
model.Lp = state.Lp;
model.dx = state.dx;
model.g = state.g;
model.feed = state.feed;
model.k_r = state.k_r;
model.C_r = state.C_r;
model.rho_r = state.rho_r;
model.tt = state.tt;
model.t_adjust = state.t_adjust;
model.t_transit = state.t_transit;
model.mCut = max([1e-4, 1e-3 * max(abs(state.Q_init), 1), 0.2]);
model.useRockAtZeroFlow = isfield(state, 'hasRockTemperatureData') && state.hasRockTemperatureData;
model.hRelax = 1 / max(state.dx, 1e-3);

xNodes = state.x(:);
sNodes = state.Lp - xNodes;
[sNodesAsc, idxNode] = sort(sNodes, 'ascend');

model.D = griddedInterpolant(sNodesAsc, state.Dp(idxNode), 'linear', 'nearest');
model.eps = griddedInterpolant(sNodesAsc, state.eps(idxNode), 'linear', 'nearest');

model.PI = griddedInterpolant(sNodesAsc, state.feedzone_PI(idxNode), 'linear', 'nearest');
model.Pres = griddedInterpolant(sNodesAsc, state.feedzone_P_res(idxNode), 'linear', 'nearest');
model.Hres = griddedInterpolant(sNodesAsc, state.feedzone_H_res(idxNode), 'linear', 'nearest');
model.Qm = griddedInterpolant(sNodesAsc, state.feedzone_Qm(idxNode), 'linear', 'nearest');
model.Hq = griddedInterpolant(sNodesAsc, state.H_q_node(idxNode), 'linear', 'nearest');

xCell = 0.5 * (state.x(1:end-1) + state.x(2:end));
sCell = state.Lp - xCell(:);
[sCellAsc, idxCell] = sort(sCell, 'ascend');
Trock = state.T_rock(:);
if isempty(Trock)
    Trock = (state.T_surface + 273.15) * ones(size(sCellAsc));
else
    Trock = Trock(idxCell);
end
model.Trock = griddedInterpolant(sCellAsc, Trock, 'linear', 'nearest');
end

function [P, H, Mdot, pBotUsed] = solveTopDown(state, model)
pBotUsed = state.P_bot;

pTop = state.P_top;
if isfinite(state.T_top_fin)
    hTop = 1e3 * IAPWS_IF97('h_pT', pTop/1e6, state.T_top_fin + 273.15);
else
    hTop = 1e3 * state.H_top;
end
mTop = max(state.Q_init, 0);

sSpan = linspace(0, state.Lp, state.n); % top -> bottom
opts = solverOptions(state);

warnState = warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
cleanupObj = onCleanup(@() warning(warnState));
[sSol, ySol] = ode15s(@(s, y) steadyRHS(s, y, model), sSpan, [pTop; hTop; mTop], opts);
clear cleanupObj

idxFinite = find(all(isfinite(ySol), 2), 1, 'last');
if isempty(idxFinite)
    error('initializeSteadyState:TopDownIntegrationFailed', ...
        'Top-down integration returned no finite states.');
end
sSol = sSol(1:idxFinite);
ySol = ySol(1:idxFinite, :);

yFull = nan(numel(sSpan), 3);
nSolved = min(numel(sSol), numel(sSpan));
yFull(1:nSolved, :) = ySol(1:nSolved, :);

if nSolved < numel(sSpan)
    yPrev = yFull(nSolved, :).';
    for k = nSolved+1:numel(sSpan)
        sPrev = sSpan(k-1);
        sNow = sSpan(k);
        ds = sNow - sPrev;

        if yPrev(3) > model.mCut
            yp = steadyRHS(sPrev, yPrev, model);
            yNow = yPrev + ds * yp;
            yNow(1) = max(yNow(1), 1e3);
            yNow(3) = max(yNow(3), 0);
            if ~isfinite(yNow(2))
                yNow(2) = yPrev(2);
            end
        else
            qMassLine = sourceMassLine(sPrev, yPrev(1), model);
            mNow = max(0, yPrev(3) - ds * qMassLine);
            base = evalLocal(sPrev, [yPrev(1); yPrev(2); 0], model);
            pPrime = base.Sgrav / max(base.A, 1e-12);
            pNow = max(1e3, yPrev(1) + ds * pPrime);
            if model.useRockAtZeroFlow
                hNow = rockEnthalpy(sNow, pNow, model);
            else
                hNow = yPrev(2);
            end
            yNow = [pNow; hNow; mNow];
        end

        if any(~isfinite(yNow))
            yNow = yPrev;
            yNow(3) = 0;
        end
        yFull(k, :) = yNow.';
        yPrev = yNow;
    end
end

if size(yFull,1) ~= numel(sSpan)
    yFull = interp1(sSol, ySol, sSpan, 'linear', 'extrap');
end

YbotToTop = yFull(end:-1:1, :).';
P = max(YbotToTop(1, :), 1e3);
H = YbotToTop(2, :);
Mdot = max(YbotToTop(3, :), 0);
if model.useRockAtZeroFlow
    sBotToTop = linspace(state.Lp, 0, state.n);
    for i = 1:state.n
        if Mdot(i) <= model.mCut
            H(i) = rockEnthalpy(sBotToTop(i), P(i), model);
        end
    end
end
end

function [P, H, Mdot, pBotUsed] = solveBottomUp(state, model)
pTopTarget = state.P_top;
pTry = state.P_bot;
maxShoot = 80;

bottomThermo = inferBottomThermoMode(state, pTry);
if bottomThermo.isTwoPhase
    fprintf('  Bottom two-phase mode detected: x=%.4f, enforcing Tsat(P_bot) coupling.\n', bottomThermo.quality);
end

sSpan = linspace(state.Lp, 0, state.n); % bottom -> top
opts = solverOptions(state);

bestY = [];
pBotUsed = pTry;

for iter = 1:maxShoot
    hBotTry = bottomEnthalpyAtPressure(state, pTry, bottomThermo);
    y0 = [pTry; hBotTry; max(state.Q_init, 0)];

    warnState = warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
    cleanupObj = onCleanup(@() warning(warnState));
    [sSol, ySol] = ode15s(@(s, y) steadyRHS(s, y, model), sSpan, y0, opts);
    clear cleanupObj

    idxFinite = find(all(isfinite(ySol), 2), 1, 'last');
    if isempty(idxFinite)
        reachedTop = false;
        pStop = pTry;
    else
        sSol = sSol(1:idxFinite);
        ySol = ySol(1:idxFinite, :);
        reachedTop = idxFinite == numel(sSpan) && ...
            sSol(end) <= sSpan(end) + 1e-10 * max(1, state.Lp);
        pStop = ySol(end, 1);
    end

    if reachedTop && isfinite(pStop) && pStop >= pTopTarget
        yFull = ySol;
        if size(yFull,1) ~= numel(sSpan)
            yFull = interp1(sSol, ySol, sSpan, 'linear', 'extrap');
        end
        if all(isfinite(yFull(:)))
            bestY = yFull.';
            pBotUsed = pTry;
            break;
        end
    end

    deficit = pTopTarget - pStop;
    if ~isfinite(deficit) || deficit <= 0
        deficit = max(1e5, 0.1 * max(pTopTarget, 1));
    end
    pTry = pTry + max(5e4, 1.2 * deficit);
end

if isempty(bestY)
    error('initializeSteadyState:BottomUpShootingFailed', ...
        'Failed to reach top boundary with P>=P_top after increasing P_bot.');
end

P = max(bestY(1, :), 1e3);
H = bestY(2, :);
Mdot = max(bestY(3, :), 0);
if model.useRockAtZeroFlow
    sBotToTop = linspace(state.Lp, 0, state.n);
    for i = 1:state.n
        if Mdot(i) <= model.mCut
            H(i) = rockEnthalpy(sBotToTop(i), P(i), model);
        end
    end
end
end

function [P, H, Mdot, pBotUsed] = solvePrescribed(state, ~)
% Prescribed-profile initial condition (IC_switch = 3).
% Pressure is taken from InitPressure.csv and temperature from
% InitTemperature.csv (both loaded into state.initProfiles in the x-position
% coordinate, 0 = bottom, Lp = top). Enthalpy follows from IAPWS h(P,T) and
% the velocity field is set so each cell carries the target mass flux Q_init.
pp = state.initProfiles.pressure;
tp = state.initProfiles.temperature;
if isempty(pp.position) || isempty(tp.position)
    error('initializeSteadyState:MissingInitProfiles', ...
        ['IC_switch = 3 requires InitPressure.csv and InitTemperature.csv ' ...
         'in %s.'], state.SimDir);
end

x = state.x(:).';
P = interpProfile(pp.position, pp.value, x);   % Pa
T = interpProfile(tp.position, tp.value, x);   % degC
P = max(P, 1e3);

H = zeros(1, state.n);
for i = 1:state.n
    H(i) = 1e3 * IAPWS_IF97('h_pT', P(i)/1e6, T(i) + 273.15);
end

Mdot = max(state.Q_init, 0) * ones(1, state.n);
pBotUsed = P(1);
end

function v = interpProfile(pos, val, xq)
% Linear interpolation of a profile onto the grid xq, sorted and de-duplicated,
% with nearest-value hold outside the sampled range.
pos = pos(:); val = val(:);
[pos, idx] = sort(pos); val = val(idx);
[pos, iu] = unique(pos, 'stable'); val = val(iu);
v = interp1(pos, val, xq(:).', 'linear', 'extrap');
v(xq < pos(1))   = val(1);
v(xq > pos(end)) = val(end);
end

function opts = solverOptions(state)
opts = odeset('RelTol', 1e-4, ...
    'AbsTol', [1e3, 1e3, 1e-3], ...
    'InitialStep', max(state.dx * 0.25, 1e-3), ...
    'MaxStep', max(state.dx, 1e-3));
end

function yp = steadyRHS(s, y, model)
y = y(:);
y(1) = max(y(1), 1e3);
y(3) = max(y(3), 0);

base = evalLocal(s, y, model);
mPrime = -base.qMassLine;

% In hydrostatic/small-flow regime the full energy equation is degenerate.
% Enforce hydrostatic pressure and relaxation to rock temperature.
if y(3) <= model.mCut
    pPrime = base.Sgrav / max(base.A, 1e-12);
    if model.useRockAtZeroFlow
        [hRock, dhdp] = rockEnthalpyAndSlope(s, y(1), model);
        hPrime = dhdp * pPrime + model.hRelax * (hRock - y(2));
    else
        hPrime = 0;
    end
    yp = [pPrime; hPrime; mPrime];
    return;
end

[dMom, dEnergy] = fluxPartials(s, y, base, model);

A11 = base.A - dMom(1);
A12 = -dMom(2);
A21 = dEnergy(1);
A22 = dEnergy(2);

b1 = dMom(3) * mPrime + base.Sfric + base.Sgrav;
b2 = base.m * model.g * base.cth + base.qHeat - base.qEnergyLine - dEnergy(3) * mPrime;

M = [A11, A12; A21, A22];
rhs = [b1; b2];

if any(~isfinite(M(:))) || any(~isfinite(rhs)) || rcond(M) < 1e-11
    pPrime = base.Sgrav / max(base.A, 1e-12);
    if y(3) <= model.mCut
        if model.useRockAtZeroFlow
            [hRock, dhdp] = rockEnthalpyAndSlope(s, y(1), model);
            hPrime = dhdp * pPrime + model.hRelax * (hRock - y(2));
        else
            hPrime = 0;
        end
    else
        hPrime = 0;
    end
else
    sol = M \ rhs;
    pPrime = sol(1);
    hPrime = sol(2);
end

if ~isfinite(pPrime)
    pPrime = base.Sgrav / max(base.A, 1e-12);
end
if ~isfinite(hPrime)
    hPrime = 0;
end

yp = [pPrime; hPrime; mPrime];
end

function [dMom, dEnergy] = fluxPartials(s, y, base, model)
dMom = zeros(1, 3);
dEnergy = zeros(1, 3);

steps = [max(100, 1e-6 * abs(y(1))), ...
         max(100, 1e-6 * abs(y(2))), ...
         max(1e-5, 1e-6 * max(1, abs(y(3))))];

for j = 1:3
    dy = steps(j);
    yP = y; yM = y;
    yP(j) = yP(j) + dy;
    yM(j) = yM(j) - dy;
    yP(1) = max(yP(1), 1e3);
    yM(1) = max(yM(1), 1e3);
    yP(3) = max(yP(3), 0);
    yM(3) = max(yM(3), 0);

    bP = evalLocal(s, yP, model);
    bM = evalLocal(s, yM, model);

    dMom(j) = (bP.momFlux - bM.momFlux) / (2 * dy);
    dEnergy(j) = (bP.energyFlux - bM.energyFlux) / (2 * dy);
end

if ~isfinite(base.momFlux)
    dMom(:) = 0;
end
if ~isfinite(base.energyFlux)
    dEnergy(:) = 0;
end
end

function base = evalLocal(s, y, model)
state = model.state;

p = max(y(1), 1e3);
h = y(2);
m = max(y(3), 0);

[D, A, rough, theta, cth] = localGeometry(s, model);
[ag, al, rg, rl, hg, hl, mu, T] = calculatePhaseProperties(p, h, state);
rho = ag * rg + al * rl;

uMix = 0;
if m > model.mCut
    uMix = m / max(A * rho, 1e-12);
end

[ug, ul] = calculatePhaseVelocities(uMix, ag, rg, rl, D, T, theta, state);

mdotG = A * ag * rg * ug;
mdotL = A * al * rl * ul;
mdotTot = mdotG + mdotL;
if abs(mdotTot) > 1e-12
    fac = m / mdotTot;
    mdotG = fac * mdotG;
    mdotL = fac * mdotL;
    ug = fac * ug;
    ul = fac * ul;
end

momFlux = mdotG * ug + mdotL * ul;
energyFlux = mdotG * (hg + 0.5 * ug^2) + mdotL * (hl + 0.5 * ul^2);

Re = max(1e-1, rho * abs(uMix) * D / max(mu, 1e-12));
fD = frictionFactor(D, Re, rough);
Sfric = A * fD * rho * uMix * abs(uMix) / (2 * max(D, 1e-12));
Sgrav = rho * model.g * A * cth;

qMassLine = sourceMassLine(s, p, model);
qEnergyLine = qMassLine * model.Hres(s);
qHeat = localHeatFlux(s, D, A, T, model);

base = struct('p', p, 'h', h, 'm', m, ...
    'D', D, 'A', A, 'cth', cth, ...
    'momFlux', momFlux, 'energyFlux', energyFlux, ...
    'Sfric', Sfric, 'Sgrav', Sgrav, ...
    'qMassLine', qMassLine, 'qEnergyLine', qEnergyLine, 'qHeat', qHeat);
end

function qMassLine = sourceMassLine(s, p, model)
if model.feed == 0
    pRes = model.Pres(s);
    if ~isfinite(pRes) || p >= pRes
        sCell = 0;
    else
        sCell = model.PI(s) * (pRes - p);
    end
    if ~isfinite(sCell) || sCell < 0
        sCell = 0;
    end
elseif model.feed == 3
    sCell = model.Qm(s);
    if sCell < 0
        sCell = 0;
    end
else
    sCell = 0;
end
qMassLine = sCell / max(model.dx, 1e-12);
end

function qHeat = localHeatFlux(s, D, A, T, model)
if model.tt < model.t_adjust + model.t_transit / 2
    qHeat = 0;
    return;
end
Trock = model.Trock(s);
H_q = model.Hq(s);
kappa = model.k_r / max(model.C_r * model.rho_r, 1e-12);
tD = kappa * model.tt / max(D^2, 1e-12);
ftd = 0.982 * log(1 + 1.81 * tD);
qHeat = A * 4 / max(D, 1e-12) * H_q * model.k_r * (T - Trock) / ...
    max(model.k_r + H_q * D * ftd / 2, 1e-12);
end

function [D, A, rough, theta, cth] = localGeometry(s, model)
D = model.D(s);
A = pi * D^2 / 4;
rough = model.eps(s);
xPos = model.Lp - s;
[theta, cth] = getGravityInclination(model.state, xPos);
end

function fD = frictionFactor(D, Re, rough)
term = (2e4 * max(rough, 0) / max(D, 1e-12) + 1e6 / max(Re, 1e-12))^(1/3);
fD = 0.0055 * (1 + term);
end

function [hRock, dhdp] = rockEnthalpyAndSlope(s, p, model)
hRock = rockEnthalpy(s, p, model);
dp = max(100, 1e-6 * abs(p));
hF = rockEnthalpy(s, p + dp, model);
hB = rockEnthalpy(s, p - dp, model);
dhdp = (hF - hB) / (2 * dp);
end

function hRock = rockEnthalpy(s, p, model)
Trock = model.Trock(s);
hRock = 1e3 * IAPWS_IF97('h_pT', max(p, 1e3)/1e6, Trock);
end

function u = velocityFromMdot(mdot, p, h, D, state, mCut)
if mdot <= mCut
    u = 0;
    return;
end
[ag, al, rg, rl] = calculatePhaseProperties(p, h, state);
rho = ag * rg + al * rl;
A = pi * D^2 / 4;
u = mdot / max(A * rho, 1e-12);
end

function mode = inferBottomThermoMode(state, pBot)
mode = struct('isTwoPhase', false, 'quality', 0.0);
TbotK = state.T_bot + 273.15;
try
    hBot = 1e3 * IAPWS_IF97('h_pT', pBot/1e6, TbotK);
    hL = 1e3 * IAPWS_IF97('hL_p', pBot/1e6);
    hV = 1e3 * IAPWS_IF97('hV_p', pBot/1e6);
catch
    return;
end

dh = max(hV - hL, 1e-6);
x = (hBot - hL) / dh;
x = min(1, max(0, x));
if x > 1e-4 && x < 1 - 1e-4
    mode.isTwoPhase = true;
    mode.quality = x;
end
end

function hBot = bottomEnthalpyAtPressure(state, pBot, mode)
if mode.isTwoPhase
    hL = 1e3 * IAPWS_IF97('hL_p', pBot/1e6);
    hV = 1e3 * IAPWS_IF97('hV_p', pBot/1e6);
    hBot = hL + mode.quality * (hV - hL);
else
    hBot = 1e3 * IAPWS_IF97('h_pT', pBot/1e6, state.T_bot + 273.15);
end
end
