function [Mliq, Mgas] = splitComponentMassByPartitionV2(Mtot, K, mLiq, mGas, gasMask)
% splitComponentMassByPartitionV2 Split total component mass by Kw = wg/wl.

[m, ~] = size(Mtot);
mLiq = mLiq(:).';
mGas = mGas(:).';
gasMask = gasMask(:);
if numel(gasMask) < m
    gasMask(end + 1:m, 1) = false;
elseif numel(gasMask) > m
    gasMask = gasMask(1:m);
end

K(~gasMask, :) = 0;
den = mLiq + K .* mGas;
den = max(den, 1e-12); % Avoid division by zero
fracL = mLiq ./ den;

Mliq = Mtot .* fracL;
Mgas = Mtot - Mliq;
Mgas(~gasMask, :) = 0;
end
