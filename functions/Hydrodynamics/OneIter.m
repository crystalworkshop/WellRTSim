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
        disp([S(:,i), Y(:,i), Y0(:,i)]);
        error('NaN/Inf in residual or Jacobian at i=%d', i);
    end
end

[dY, id] = TDMA(A, C, B, S, n);
if id == 1
    disp('TDMA failed, reducing timestep');
    Y = Y0;
    tol1 = Inf;
    tol2 = Inf;
    state.tt = state.tt - state.dt;
    state.dt = state.dt / 2;
    return
end

if min(state.Dp) < 0.01
    warning('pipe diameter is too small');
    state.cancelFlag = true;
    Y = Y0;
    tol1 = Inf;
    tol2 = Inf;
    return;
end

Y = Y + dY;
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

sc = [max(1e3, 1e-6*abs(Zt(1))); ...
      max(1e3, 1e-6*abs(Zt(2))); ...
      max(1e-4, 1e-3*abs(Zt(3)))];

for j = 1:k
    Zp = Zt;
    h = sc(j);
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
