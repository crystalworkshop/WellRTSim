function CliqNew = chemistryTransportStepV2(state, Cliq, geom)
% chemistryTransportStepV2 Transport liquid-phase mass fractions implicitly.

[m, n] = size(Cliq);
dt = state.dt;

K = state.chem.partitionCoefficients;
gasMask = state.chem.gasMask(:);
K(~gasMask, :) = 0;

mLiq = geom.mLiquid;
mGas = geom.mGas;
qLiqFace = state.Q_l_face(1:n + 1);
qGasFace = state.Q_v_face(1:n + 1);
wLiqFace = state.w_l_face(1:n + 1);
wGasFace = state.w_v_face(1:n + 1);

dmSource = zeros(1, n);
if state.chem_source ~= 0
    sourceLiquidFrac = ones(1, n);
    activeSources = state.feedzone_cells(:).' == 1 & (state.Q_mass(1:n) > 0);
    if any(activeSources)
        idx = find(activeSources);
        [alphaG, alphaL, rhoG, rhoL] = calculatePhaseProperties( ...
            state.feedzone_P_res(idx), state.feedzone_H_res(idx), state);
        qGas = alphaG(:).' .* rhoG(:).' ./ (alphaG(:).' .* rhoG(:).' + alphaL(:).' .* rhoL(:).');
        sourceLiquidFrac(idx) = 1 - qGas;
    end
    dmSource = state.Q_mass(1:n) .* sourceLiquidFrac .* (geom.Acell * state.dx) * dt;
end
sourceMass = state.inlet_massfrac(:) * dmSource;
liqBoundary = state.inlet_massfrac(:);

accum = mLiq + K .* mGas;
rhs = accum .* Cliq + sourceMass;

leftCoeffSelf = zeros(m, n);
leftCoeffSelf(:, 1) = qLiqFace(1) * (1 - wLiqFace(1)) + qGasFace(1) * (1 - wGasFace(1)) .* K(:, 1);
if n > 1
    leftCoeffSelf(:, 2:n) = ...
        qLiqFace(2:n) .* (1 - wLiqFace(2:n)) + ...
        qGasFace(2:n) .* (1 - wGasFace(2:n)) .* K(:, 2:n);
end

rightCoeffSelf = zeros(m, n);
if n > 1
    rightCoeffSelf(:, 1:n-1) = ...
        qLiqFace(2:n) .* wLiqFace(2:n) + ...
        qGasFace(2:n) .* wGasFace(2:n) .* K(:, 1:n-1);
end
rightCoeffSelf(:, n) = ...
    qLiqFace(n + 1) * wLiqFace(n + 1) + ...
    qGasFace(n + 1) * wGasFace(n + 1) .* K(:, n);

a = zeros(m, n);
if n > 1
    a(:, 2:n) = -dt * ( ...
        qLiqFace(2:n) .* wLiqFace(2:n) + ...
        qGasFace(2:n) .* wGasFace(2:n) .* K(:, 1:n-1));
end

b = zeros(m, n);
if n > 1
    b(:, 1:n-1) = dt * ( ...
        qLiqFace(2:n) .* (1 - wLiqFace(2:n)) + ...
        qGasFace(2:n) .* (1 - wGasFace(2:n)) .* K(:, 2:n));
end

c = accum + dt * (rightCoeffSelf - leftCoeffSelf);
rhs(:, 1) = rhs(:, 1) + dt * ( ...
    qLiqFace(1) * wLiqFace(1) * liqBoundary);
rhs(:, n) = rhs(:, n) - dt * ( ...
    qLiqFace(n + 1) * (1 - wLiqFace(n + 1)) * liqBoundary);

CliqNew = tridiag_solve_batch(a, c, b, rhs);
end

function x = tridiag_solve_batch(a, c, b, r)
[m, n] = size(c);
x = zeros(m, n);
cp = zeros(m, n);
dp = zeros(m, n);

cp(:, 1) = b(:, 1) ./ c(:, 1);
dp(:, 1) = r(:, 1) ./ c(:, 1);

for i = 2:n
    den = c(:, i) - a(:, i) .* cp(:, i - 1);
    cp(:, i) = b(:, i) ./ den;
    dp(:, i) = (r(:, i) - a(:, i) .* dp(:, i - 1)) ./ den;
end

x(:, n) = dp(:, n);
for i = n-1:-1:1
    x(:, i) = dp(:, i) - cp(:, i) .* x(:, i + 1);
end
end
