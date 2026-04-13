function [Cliq, Mliq, mWater] = phaseMolesToLiquidMassfracV2(Nliq, mLiq, molarMassKg)
% phaseMolesToLiquidMassfracV2 Convert liquid-phase moles to mass fractions.

mLiq = mLiq(:).';
molarMassKg = molarMassKg(:);

Mliq = Nliq .* molarMassKg;
mWater = mLiq - sum(Mliq, 1);
if any(mWater <= 0)
    error('phaseMolesToLiquidMassfracV2:InvalidWaterMass', ...
        'Liquid component masses exceed the liquid-phase mass.');
end
Cliq = Mliq ./ mLiq;
end
