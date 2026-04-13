function balance = computeCalciumMassBalanceV2(state, geom, NliqBefore, NgasBefore, kinMolesBefore, NliqAdv, NgasAdv)
% computeCalciumMassBalanceV2 Compute per-step Ca mass balance diagnostics.
% Hydro/domain terms are over the hydrodynamic timestep.
% Chem terms are over the reaction/flight step used in PHREEQC.

caIdx = state.chem.caComponentIndex;
mwCa = state.molar_mass(caIdx);
n = geom.n;
dt = state.dt;
liqBoundaryCa = state.inlet_massfrac(caIdx);
cellCa = state.C(caIdx, 1:n);

dmSource = compute_liquid_source_mass_step(state, geom);
caSourceSignedKg = liqBoundaryCa .* dmSource;
sourceInKg = sum(max(caSourceSignedKg, 0));
sourceOutKg = sum(max(-caSourceSignedKg, 0));

qLiqFace = state.Q_l_face(1:n + 1);
wLiqFace = state.w_l_face(1:n + 1);
boundaryInKg = ...
    max(dt * qLiqFace(1) * wLiqFace(1) * liqBoundaryCa, 0) + ...
    max(-dt * qLiqFace(n + 1) * (1 - wLiqFace(n + 1)) * liqBoundaryCa, 0);
boundaryOutKg = ...
    max(-dt * qLiqFace(1) * (1 - wLiqFace(1)) * cellCa(1), 0) + ...
    max(dt * qLiqFace(n + 1) * wLiqFace(n + 1) * cellCa(n), 0);

preFluidKg = compute_component_mass_kg(NliqBefore, NgasBefore, mwCa, caIdx);
postFluidKg = compute_component_mass_kg(state.chem.N_liq, state.chem.N_gas, mwCa, caIdx);
chemInKg = compute_component_mass_kg(NliqAdv, NgasAdv, mwCa, caIdx);

deltaKin = state.chem.kinMoles - kinMolesBefore;
caStoich = state.chem.mineralCalciumStoich(:);
solidDeltaMolCaHydro = sum(reshape(deltaKin .* caStoich, 1, []));
reactionScale = reshape(state.chem.lastReactionScale, 1, []);
deltaKinReaction = deltaKin ./ repmat(reactionScale, size(deltaKin, 1), 1);
solidDeltaMolCaReaction = sum(reshape(deltaKinReaction .* caStoich, 1, []));

balance = struct( ...
    'componentIndex', caIdx, ...
    'inKg', sourceInKg + boundaryInKg, ...
    'outKg', sourceOutKg + boundaryOutKg, ...
    'precipKg', solidDeltaMolCaReaction * mwCa, ...
    'precipHydroKg', solidDeltaMolCaHydro * mwCa, ...
    'fluidDeltaKg', postFluidKg - preFluidKg, ...
    'residualKg', sourceInKg + boundaryInKg - sourceOutKg - boundaryOutKg - solidDeltaMolCaHydro * mwCa - (postFluidKg - preFluidKg), ...
    'chemInKg', chemInKg, ...
    'chemOutKg', postFluidKg, ...
    'chemResidualKg', chemInKg - postFluidKg - solidDeltaMolCaReaction * mwCa, ...
    'sourceInKg', sourceInKg, ...
    'sourceOutKg', sourceOutKg, ...
    'boundaryInKg', boundaryInKg, ...
    'boundaryOutKg', boundaryOutKg);
end

function totalKg = compute_component_mass_kg(Nliq, Ngas, molarMassKg, idx)
totalKg = sum(Nliq(idx, :) + Ngas(idx, :)) * molarMassKg;
end

function dmSource = compute_liquid_source_mass_step(state, geom)
n = geom.n;
dmSource = zeros(1, n);
if state.chem_source == 0
    return
end

sourceLiquidFrac = ones(1, n);
activeSources = state.feedzone_cells(:).' == 1 & (state.Q_mass(1:n) > 0);
if any(activeSources)
    idx = find(activeSources);
    [alphaG, alphaL, rhoG, rhoL] = calculatePhaseProperties( ...
        state.feedzone_P_res(idx), state.feedzone_H_res(idx), state);
    qGas = alphaG(:).' .* rhoG(:).' ./ (alphaG(:).' .* rhoG(:).' + alphaL(:).' .* rhoL(:).');
    sourceLiquidFrac(idx) = 1 - qGas;
end

dmSource = state.Q_mass(1:n) .* sourceLiquidFrac .* (geom.Acell * state.dx) * state.dt;
end
