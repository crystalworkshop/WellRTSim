function T = computeTemperatureProfile(Y, state)
% computeTemperatureProfile Compute temperature profile [K] from [P; H; U].

P = Y(1, :);
H = Y(2, :);

if state.interp
    T = state.Temp(H, P);
else
    T = IAPWS_IF97('T_ph', P / 1e6, H / 1e3);
end
end
