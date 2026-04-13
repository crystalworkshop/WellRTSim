function state = chemistryStepV2(state)
% chemistryStepV2 One chemistry update with optional semi-implicit coupling:
%   hydrodynamic fields -> transport -> PHREEQC closure, iterated within dt.

if ~state.chem.initialized
    state = chemistryInitializeStateV2(state);
end

geom = chemistryPhaseGeometryV2(state);
NliqBefore = state.chem.N_liq;
NgasBefore = state.chem.N_gas;
kinMolesBefore = state.chem.kinMoles;
liqResidence = geom.mLiquid ./ max(abs(state.Q_l(1:geom.n)), eps);
reactionDt = min(state.dt * ones(1, geom.n), liqResidence);
reactionScale = state.dt ./ max(reactionDt, eps);

% Transport uses K = w_g / w_l in mass fractions.
% Mode 2 refreshes K analytically before transport.
% Mode 1 keeps the K obtained from the previous PHREEQC gas-phase split.
if state.chem.useAnalyticalGasPartition
    K = partitionCoefficients(state, geom.T_K, geom.rho_l, geom.rho_g, ...
        geom.alpha_g, geom.alpha_l);
    K(~state.chem.gasMask(:), :) = 0;
    state.chem.partitionCoefficients = K;
end

if should_use_semi_implicit_coupling(state)
    [Nliq, Ngas, out, semi, NliqAdv, NgasAdv] = run_semi_implicit_coupling(state, geom, reactionDt, reactionScale);
    state.chem.lastSemiImplicitIterations = semi.iterations;
    state.chem.lastSemiImplicitResidualK = semi.residualK;
    state.chem.lastSemiImplicitResidualMoles = semi.residualMoles;
    state.chem.lastSemiImplicitConverged = semi.converged;
else
    [NliqAdv, NgasAdv] = prepare_advected_phase_moles( ...
        state, geom, state.chem.partitionCoefficients);
    [Nliq, Ngas, out] = runPhreeqcStepV2(state, NliqAdv, NgasAdv, reactionDt, reactionScale, true, geom);
    state.chem.lastSemiImplicitIterations = 0;
    state.chem.lastSemiImplicitResidualK = NaN;
    state.chem.lastSemiImplicitResidualMoles = NaN;
    state.chem.lastSemiImplicitConverged = false;
end

state.chem.lastReactionTimeS = out.reactionDtS(:).';
state.chem.lastReactionScale = out.reactionScale(:).';
state.chem.lastPhreeqcScript = out.script;
state = updateChemistryDerivedV2( ...
    state, geom, Nliq, Ngas, ...
    out.pH, out.si, out.kinMoles, out.partitionCoefficients);
state.chem.lastCalciumBalance = computeCalciumMassBalanceV2( ...
    state, geom, NliqBefore, NgasBefore, kinMolesBefore, NliqAdv, NgasAdv);
end

function tf = should_use_semi_implicit_coupling(state)
tf = state.chem.semiImplicitEnabled ...
    && ~state.chem.useAnalyticalGasPartition ...
    && any(state.chem.gasMask(:));
end

function [Nliq, Ngas, out, info, NliqAdvFinal, NgasAdvFinal] = run_semi_implicit_coupling(state, geom, reactionDt, reactionScale)
maxIter = state.chem.semiImplicitMaxIter;
relax = state.chem.semiImplicitRelaxation;
tolK = state.chem.semiImplicitKTolerance;
tolMoles = state.chem.semiImplicitMolesTolerance;

Kiter = sanitize_partition_coefficients(state.chem.partitionCoefficients, state.chem.gasMask);
NliqPrev = [];
NgasPrev = [];

Nliq = state.chem.N_liq;
Ngas = state.chem.N_gas;
out = struct();
resK = NaN;
resMoles = NaN;
converged = false;
NliqAdvFinal = state.chem.N_liq;
NgasAdvFinal = state.chem.N_gas;

for iter = 1:maxIter
    [NliqAdv, NgasAdv] = prepare_advected_phase_moles(state, geom, Kiter);
    NliqAdvFinal = NliqAdv;
    NgasAdvFinal = NgasAdv;

    stateIter = state;
    stateIter.chem.partitionCoefficients = Kiter;

    [NliqRaw, NgasRaw, outRaw] = runPhreeqcStepV2( ...
        stateIter, NliqAdv, NgasAdv, reactionDt, reactionScale, true, geom);
    Kraw = sanitize_partition_coefficients(outRaw.partitionCoefficients, state.chem.gasMask);

    if iter > 1
        resK = max_relative_change(Kiter(state.chem.gasMask(:), :), ...
            Kraw(state.chem.gasMask(:), :));
        resMoles = max_relative_change(NliqPrev + NgasPrev, NliqRaw + NgasRaw);
        if resK <= tolK && resMoles <= tolMoles
            Nliq = NliqRaw;
            Ngas = NgasRaw;
            out = outRaw;
            converged = true;
            break
        end
    end

    Nliq = NliqRaw;
    Ngas = NgasRaw;
    out = outRaw;
    NliqPrev = NliqRaw;
    NgasPrev = NgasRaw;

    if iter < maxIter
        Kiter = relax_partition_coefficients(Kiter, Kraw, relax, state.chem.gasMask);
    end
end

info = struct();
info.iterations = iter;
info.residualK = resK;
info.residualMoles = resMoles;
info.converged = converged || maxIter == 1;
end

function [NliqAdv, NgasAdv] = prepare_advected_phase_moles(state, geom, Ktransport)
Ktransport = sanitize_partition_coefficients(Ktransport, state.chem.gasMask);

stateTransport = state;
stateTransport.chem.partitionCoefficients = Ktransport;
CliqAdv = chemistryTransportStepV2(stateTransport, state.C, geom);
CliqAdv = max(CliqAdv, 0);

[NliqAdv, ~, ~] = liquidMassfracToMolesV2(CliqAdv, geom.mLiquid, state.molar_mass);
MgasAdv = Ktransport .* CliqAdv .* geom.mGas;
NgasAdv = MgasAdv ./ state.molar_mass(:);
NgasAdv(~isfinite(NgasAdv) | NgasAdv < 0) = 0;
end

function K = sanitize_partition_coefficients(K, gasMask)
K(~gasMask(:), :) = 0;
end

function Knew = relax_partition_coefficients(Kold, Kraw, relax, gasMask)
Knew = sanitize_partition_coefficients(Kold, gasMask);
Kraw = sanitize_partition_coefficients(Kraw, gasMask);
mask = repmat(gasMask(:), 1, size(Knew, 2));
Knew(mask) = (1 - relax) .* Knew(mask) + relax .* Kraw(mask);
Knew(~mask) = 0;
end

function r = max_relative_change(oldVal, newVal)
scale = max(max(abs(oldVal), abs(newVal)), 1e-12);
r = max(abs(newVal - oldVal) ./ scale, [], 'all');
end
