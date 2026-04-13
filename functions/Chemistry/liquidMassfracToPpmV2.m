function ppm = liquidMassfracToPpmV2(Cliq)
% liquidMassfracToPpmV2 Convert liquid-phase mass fractions to ppm (mg/kgw).

waterMassfrac = 1 - sum(Cliq, 1);
ppm = 1e6 * Cliq ./ waterMassfrac;
end
