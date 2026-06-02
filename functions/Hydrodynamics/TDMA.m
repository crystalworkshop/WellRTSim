function [dy, id] = TDMA(ax, cx, bx, fx, n)
k = size(fx,1);
Psi = zeros(k, k, n);
teta = zeros(k, n);
dy = zeros(k, n);
id = 0;

[M, ok] = safeSolveMatrix(cx(:,:,1), 1);
if ~ok, id = 1; return; end
Psi(:,:,1) = M * bx(:,:,1);
teta(:,1) = M * fx(:,1);

for i = 1:n-1
    [M, ok] = safeSolveMatrix(cx(:,:,i) - ax(:,:,i)*Psi(:,:,i), i);
    if ~ok, id = 1; return; end
    Psi(:,:,i+1) = M * bx(:,:,i);
    teta(:,i+1) = M * (ax(:,:,i)*teta(:,i) + fx(:,i));
end

[M, ok] = safeSolveMatrix(cx(:,:,n) - ax(:,:,n)*Psi(:,:,n), n);
if ~ok, id = 1; return; end
dy(:,n) = M * (fx(:,n) + ax(:,:,n)*teta(:,n));

for i = n-1:-1:1
    dy(:,i) = Psi(:,:,i+1)*dy(:,i+1) + teta(:,i+1);
end
end

function [Minv, ok] = safeSolveMatrix(A, idx)
ok = true;
if rcond(A) < 1e-14 || any(~isfinite(A), 'all')
    % Ill-conditioned block: signal failure quietly; OneStep retries with a
    % smaller dt (rejected iterates can produce singular blocks).
    Minv = zeros(size(A));
    ok = false;
    return;
end
Minv = A \ eye(size(A));
end
