function state = updateChemistryDerivedV2(state, geom, Nliq, Ngas, pH, siKin, kinMoles, Kchem)
% updateChemistryDerivedV2 Refresh state fields derived from phase moles.

[Cliq, ~, ~] = phaseMolesToLiquidMassfracV2(Nliq, geom.mLiquid, state.molar_mass);
Mgas = Ngas .* state.molar_mass(:);
if nargin >= 8
    K = Kchem;
else
    K = zeros(size(Cliq));
    gasMask = state.chem.gasMask(:);
    gasMaskMat = repmat(gasMask, 1, size(Cliq, 2));
    mGasMat = repmat(geom.mGas, size(Cliq, 1), 1);
    validK = gasMaskMat & (Cliq > 0) & (mGasMat > 0);
    K(validK) = Mgas(validK) ./ (mGasMat(validK) .* Cliq(validK));
end

state.C = Cliq;
state.chem.N_liq = Nliq;
state.chem.N_gas = Ngas;
state.chem.partitionCoefficients = K;
state.chem.pH = pH(:).';

kinNames = string(state.chem.mineralNames(:));
nKin = numel(kinNames);
n = state.n;

state.chem.saturationIndices = siKin;
state.chem.kinMoles = kinMoles;
state.chem.scaleMineralNames = kinNames;
state.chem.scaleMolarMass = state.chem.mineralMolarMass(:);
state.chem.scaleDensity = state.chem.mineralDensity(:);

D0 = state.Dp0(:).';
Dprev = state.Dp(:).';
prevTotalThickness = reshape(state.chem.totalScaleThickness, 1, []);

areaWall = pi * D0 * state.dx;
validArea = areaWall > 0;

totalThickness = zeros(1, n);
thickness = zeros(nKin, n);
arealMoles = zeros(nKin, n);

for i = 1:nKin
    mmKg = state.chem.mineralMolarMass(i);
    rho = state.chem.mineralDensity(i);
    if rho <= 0
        continue;
    end
    vol = state.chem.kinMoles(i, :) * mmKg / rho;
    validVol = validArea & vol >= 0;
    thickness(i, validVol) = vol(validVol) ./ areaWall(validVol);
    validMoles = validArea & state.chem.kinMoles(i, :) >= 0;
    arealMoles(i, validMoles) = state.chem.kinMoles(i, validMoles) ./ areaWall(validMoles);
    totalThickness = totalThickness + thickness(i, :);
end

state.chem.mineralThickness = thickness;
state.chem.mineralArealMoles = arealMoles;
state.chem.totalScaleThickness = totalThickness(:);
state.delta = totalThickness(:);

dThickness = totalThickness - prevTotalThickness;
Dnew = Dprev - 2 * dThickness;
Dnew = min(D0, Dnew);
Dnew = max(1e-6, Dnew);
state.Dp = Dnew(:);
state = updateWellRoughness(state);
end
