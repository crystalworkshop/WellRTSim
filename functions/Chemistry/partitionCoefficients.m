function K = partitionCoefficients(state, T_K, rhoW, rhoV, alpha_g, alpha_l)
% partitionCoefficients Legacy analytical gas partition model.
% Returns K = w_g/w_l for each component at each cell.

m = numel(state.chemNames);
n = numel(T_K);

K = zeros(m, n);
[KcCO2, KcH2S] = KD_CO2_H2S(T_K, rhoW, rhoV, alpha_g, alpha_l);
phaseMask = alpha_g > 0;

co2Idx = findComponentIndex(state.chemNames, ...
    {"CO2", "C", "C(4)", "Carbon Dioxide", "CarbonDioxide", "Carbon"});
h2sIdx = findComponentIndex(state.chemNames, ...
    {"H2S", "S(-2)", "Hydrogen Sulfide", "HydrogenSulfide", "Sulfide", "Sulphide"});

if co2Idx > 0
    K(co2Idx, phaseMask) = KcCO2(phaseMask);% .* rhoW(phaseMask) ./ rhoV(phaseMask);
end
if h2sIdx > 0
    K(h2sIdx, phaseMask) = KcH2S(phaseMask);% .* rhoW(phaseMask) ./ rhoV(phaseMask);
end
end

function idx = findComponentIndex(names, aliases)
idx = 0;
nameTokens = upper(string(names(:)));
aliasTokens = upper(string(aliases(:)));
aliasTokens = aliasTokens(aliasTokens ~= "");

for j = 1:numel(aliasTokens)
    match = find(nameTokens == aliasTokens(j), 1);
    if ~isempty(match)
        idx = match;
        return
    end
end
end
