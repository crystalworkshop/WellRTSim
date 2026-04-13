function [Nliq, Mliq, mWater] = liquidMassfracToMolesV2(Cliq, mLiq, molarMassKg)
% liquidMassfracToMolesV2 Convert liquid-phase mass fractions to moles.

mLiq = mLiq(:).';
molarMassKg = molarMassKg(:);

Mliq = Cliq .* mLiq;
mWater = mLiq - sum(Mliq, 1);
if any(mWater <= 0)
    error('liquidMassfracToMolesV2:InvalidWaterMass', ...
        'Liquid mass fractions leave no water mass for molality conversion.');
end
Nliq = Mliq ./ molarMassKg;
end
