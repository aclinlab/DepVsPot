function etaStar = etaStarNewton(xi, xj, yj, tol, maxIter)
% Exact eta* for f(eta) = sum xi.*exp(-eta*xi) + sum xj.*(exp(-eta*yj) + exp(-eta*xj) - 1).
% Newton from eta = 0: f is convex and decreasing, so iterates rise monotonically to the root.
if nargin < 4, tol = 1e-12; end
if nargin < 5, maxIter = 100; end
z  = [xi(:); xj(:)];
xj = xj(:);  yj = yj(:);
Sx = sum(xj);
if Sx == 0
    etaStar = Inf;      % no overlapping pairs: f > 0 for every finite eta
    return
end
eta = 0;
for k = 1:maxIter
    ez = exp(-eta*z);
    ey = exp(-eta*yj);
    f  = sum(z.*ez) + sum(xj.*ey) - Sx;
    fp = -sum(z.^2.*ez) - sum(xj.*yj.*ey);
    step = -f/fp;
    eta = eta + step;
    if abs(step) <= tol*max(eta, 1)
        break
    end
end
etaStar = eta;
end