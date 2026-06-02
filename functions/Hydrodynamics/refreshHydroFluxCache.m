function state = refreshHydroFluxCache(state, Y)
% refreshHydroFluxCache Re-evaluate RHS_v2 on the converged state to cache the
% per-face phase mass fluxes (Q_l_face/Q_v_face) and the per-phase upwind
% weights (w_l_face/w_v_face) consumed by the chemistry transport.

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
