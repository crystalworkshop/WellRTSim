function [theta, cth] = getGravityInclination(state, position)
% Return local deviation angle from vertical and the signed gravity projection.
% position is the simulation coordinate along the well, with 0 at the bottom.

cth = state.wellInterpolants.gravityProjectionX(position);
cth = max(-1, min(1, cth));
theta = acos(cth);
end
