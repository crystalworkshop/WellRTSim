function [Y, tol1, tol2, state] = OneIter(Y, Y0, state)
[k, n] = size(Y);
A = zeros(k, k, n);
B = zeros(k, k, n);
C = zeros(k, k, n);
S = zeros(k, n);

for i = 1:n
    iL = max(1, i-1);
    iR = min(n, i+1);
    Z1 = Y(:, iL);
    Z2 = Y(:, i);
    Z3 = Y(:, iR);

    [S(:,i), state] = RHS_v2(Z3, Z2, Z1, Y0, state, i);

    if i > 1
        A(:,:,i) = Jacob(Z3, Z2, Z1, S(:,i), Y0, state, i, 1);
    end
    C(:,:,i) = -Jacob(Z3, Z2, Z1, S(:,i), Y0, state, i, 2);
    if i < n
        B(:,:,i) = Jacob(Z3, Z2, Z1, S(:,i), Y0, state, i, 3);
    end

    if any(~isfinite(S(:,i))) || any(~isfinite(A(:,:,i)),'all') || ...
       any(~isfinite(B(:,:,i)),'all') || any(~isfinite(C(:,:,i)),'all')
        % Non-finite Newton iterate: signal failure so OneStep rejects the
        % step and retries with a smaller dt (no hard error / no spam).
        Y = Y0; tol1 = Inf; tol2 = Inf;
        return;
    end
end

% Equilibrate the Newton system before the solve. The four residuals and
% four primaries span ~6 orders of magnitude (energy ~1e5, mass/slip ~1;
% P,H ~1e6, FVl,uv ~1), so the raw block matrix has cond ~1e12 purely from
% units and rcond sits on the TDMA singular cutoff. Row+column scaling
% (D_r * M * D_c) recovers ~6 orders of conditioning; we solve for the
% scaled correction and unscale: dY = D_c * y.
cS = [1e6; 1e6; 1; 1];        % column scales: P[Pa] H[J/kg] FVl[m/s] uv[m/s]
S_unscaled = S;
for i = 1:n
    rowmax = max([abs(A(:,:,i)), abs(C(:,:,i)), abs(B(:,:,i))], [], 2);
    rowmax(rowmax == 0) = 1;
    rS = 1 ./ rowmax;
    A(:,:,i) = (rS .* A(:,:,i)) .* cS';
    C(:,:,i) = (rS .* C(:,:,i)) .* cS';
    B(:,:,i) = (rS .* B(:,:,i)) .* cS';
    S(:,i)   = rS .* S(:,i);
end

[dY, id] = TDMA(A, C, B, S, n);
if id == 1
    disp('TDMA failed, reducing timestep');
    Y = Y0;
    tol1 = Inf;
    tol2 = Inf;
    state.dt=state.dt / state.dt_increment; % back off timestep for retry
    return
end
dY = cS .* dY;                % unscale the correction
S = S_unscaled;               % report tol2 on the physical residual

if min(state.Dp) < 0.01
    warning('pipe diameter is too small');
    state.cancelFlag = true;
    Y = Y0;
    tol1 = Inf;
    tol2 = Inf;    
    return;
end

Y = Y + dY;
Y(1, :) = max(Y(1, :), state.P_atm);   % pressure cannot drop below atmospheric
nz = (Y ~= 0);
rel = dY.^2;
rel(nz) = (dY(nz)./Y(nz)).^2;
tol1 = sqrt(sum(rel(:)));
tol2 = sqrt(sum(S(:).^2));

if state.cans_dial.CancelRequested
    state.cancelFlag = true;
    return;
end

if state.tt < state.t_adjust && n > 1
    state.Q_mass(1) = state.Q_mass(2);
    state.Q_l(1) = state.Q_l(2);
    state.Q_v(1) = state.Q_v(2);
end
end

function Jac = Jacob(Z3, Z2, Z1, S, Y0, state, i, which)
k = numel(Z2);
Jac = zeros(k, k);

switch which
    case 1
        Zt = Z1;
    case 2
        Zt = Z2;
    otherwise
        Zt = Z3;
end

% Perturbation scales per primary: P [Pa], H [J/kg], FVl [m/s], uv [m/s].
sc = [max(1e3, 1e-6*abs(Zt(1))); ...
      max(1e3, 1e-6*abs(Zt(2))); ...
      max(1e-4, 1e-3*abs(Zt(3))); ...
      max(1e-4, 1e-3*abs(Zt(4)))];

% Phase-aware perturbation direction (Tonkin 2023, Section 4.6): step P and H
% AWAY from the phase boundary so the finite-difference gradient does not cross
% saturation (which gives spurious dalpha_g/dH and corrupts the Jacobian at
% flashing cells). Velocity primaries have no phase boundary -> step forward.
Sv_t = calculatePhaseProperties(Zt(1), Zt(2), state);
sgn = ones(k, 1);
if Sv_t <= 1e-6              % single-phase liquid: raise P, lower H (stay liquid)
    sgn(1) = +1; sgn(2) = -1;
elseif Sv_t >= 1 - 1e-6      % single-phase vapour: lower P, raise H (stay vapour)
    sgn(1) = -1; sgn(2) = +1;
else                          % two-phase: step toward the interior of the dome
    sgn(1) = +1;
    if Sv_t < 0.5, sgn(2) = +1; else, sgn(2) = -1; end
end

for j = 1:k
    Zp = Zt;
    h = sgn(j) * sc(j);
    Zp(j) = Zt(j) + h;
    switch which
        case 1
            Sp = RHS_v2(Z3, Z2, Zp, Y0, state, i);
        case 2
            Sp = RHS_v2(Z3, Zp, Z1, Y0, state, i);
        otherwise
            Sp = RHS_v2(Zp, Z2, Z1, Y0, state, i);
    end
    Jac(:,j) = (Sp - S) / h;
end
end
