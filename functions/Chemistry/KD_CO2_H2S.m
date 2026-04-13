function [Kc_CO2, Kc_H2S] = KD_CO2_H2S(TK, rhoL, ~, ~, ~)
% KD_CO2_H2S Legacy Fernandez-Prini / Alvarez gas partition correlation.
% Returns Kc = Cv/Cl for CO2 and H2S.

Tc1 = 647.096;
rhoc = 322.0;

tau = 1 - TK ./ Tc1;
f = rhoL ./ rhoc - 1;
q = -0.023767;

ECO2 = 1672.9376;  FCO2 = 28.1751;  GCO2 = -112.4619; HCO2 = 85.3807;
EH2S = 1319.1205;  FH2S = 14.1571;  GH2S = -46.8361;  HH2S = 33.2266;

ex = exp((273.15 - TK) ./ 100);

lnKD_CO2 = q .* FCO2 + (ECO2 ./ TK) .* f + ...
    (FCO2 + GCO2 .* tau .^ (2 / 3) + HCO2 .* tau) .* ex;
lnKD_H2S = q .* FH2S + (EH2S ./ TK) .* f + ...
    (FH2S + GH2S .* tau .^ (2 / 3) + HH2S .* tau) .* ex;

KD_CO2 = exp(lnKD_CO2);
KD_H2S = exp(lnKD_H2S);

Kc_CO2 = KD_CO2;
Kc_H2S = KD_H2S;
end
