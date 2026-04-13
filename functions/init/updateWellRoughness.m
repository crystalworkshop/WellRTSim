function state = updateWellRoughness(state)
% Build the effective wall roughness profile from base casing roughness and scale thickness.

if isfield(state, 'eps_base') && ~isempty(state.eps_base)
    rough = state.eps_base(:);
elseif isfield(state, 'eps') && ~isempty(state.eps)
    rough = state.eps(:);
else
    state.eps = [];
    return;
end

epsScale = 0;
if isfield(state, 'eps_scale') && ~isempty(state.eps_scale) && isfinite(state.eps_scale)
    epsScale = state.eps_scale;
end
if epsScale <= 0
    state.eps = rough;
    return;
end

scaleThickness = zeros(size(rough));

if isfield(state, 'scaleProfile') && isstruct(state.scaleProfile) && ...
        isfield(state.scaleProfile, 'position') && isfield(state.scaleProfile, 'diameter') && ...
        ~isempty(state.scaleProfile.position) && ~isempty(state.scaleProfile.diameter) && ...
        isfield(state, 'Dp0') && ~isempty(state.Dp0)
    pos = state.scaleProfile.position(:);
    dia = state.scaleProfile.diameter(:);
    valid = isfinite(pos) & isfinite(dia);
    if any(valid)
        diaGrid = interp1(pos(valid), dia(valid), state.x(:), 'nearest', nan);
        baseDia = state.Dp0(:);
        thickness0 = 0.5 * (baseDia - diaGrid);
        thickness0(~isfinite(thickness0)) = 0;
        thickness0 = max(thickness0, 0);
        scaleThickness = max(scaleThickness, thickness0);
    end
end

if isfield(state, 'chem') && isstruct(state.chem) && ...
        isfield(state.chem, 'totalScaleThickness') && ~isempty(state.chem.totalScaleThickness)
    totalThickness = state.chem.totalScaleThickness(:);
    if numel(totalThickness) == numel(rough)
        totalThickness(~isfinite(totalThickness)) = 0;
        totalThickness = max(totalThickness, 0);
        scaleThickness = max(scaleThickness, totalThickness);
    end
end

state.eps = max(rough, min(scaleThickness, epsScale));
end
