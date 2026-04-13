function state = chemistryInitializeStateV2(state)
% chemistryInitializeStateV2 Build initial split-phase chemical inventory.

geom = chemistryPhaseGeometryV2(state);
n = state.n;
m = state.mChem;

Cliq0 = repmat(state.massfrac_init(:), 1, n);
[Nliq0, ~, ~] = liquidMassfracToMolesV2(Cliq0, geom.mLiquid, state.molar_mass);
Ngas0 = zeros(m, n);

state = ensureIPhreeqcHandleV2(state);
% [Nliq, Ngas, out] = runPhreeqcStepV2(state, Nliq0, Ngas0, 0, 1, false, geom);
% state.chem.lastPhreeqcScript = out.script;
% state = updateChemistryDerivedV2( ...
%     state, geom, Nliq, Ngas, ...
%     out.pH, out.si, out.kinMoles, out.partitionCoefficients);

state.chem.initialized = true;
state.massfrac_init = state.chem.initialMassfrac(:).';
state.chem.elementInitialMass = state.C;
end
