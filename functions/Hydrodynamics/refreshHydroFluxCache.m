function state = refreshHydroFluxCache(state, Y)
% refreshHydroFluxCache Re-evaluate RHS_v2 on the converged state to cache phase fluxes.

n = state.n;
state.Q_v_face(:) = 0;
state.Q_l_face(:) = 0;
state.w_v_face(:) = 0;
state.w_l_face(:) = 0;

for i = 1:n
    iL = max(1, i - 1);
    iR = min(n, i + 1);
    [~, state] = RHS_v2(Y(:, iR), Y(:, i), Y(:, iL), Y, state, i);
end
end
