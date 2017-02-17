% Copyright 2013 Clemens Kirisits and Lukas Lang
%
% This file is part of OFD.
%
%    OFD is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFD is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFD.  If not, see <http://www.gnu.org/licenses/>.
function [c, Y] = surffit(N, X, alpha, s)
%SURFFIT Fits a sphere-like surface to data.
%
%   [c, Y] = SURFFIT(N, X, alpha, s) takes data X in R^3 and fits a 
%   sphere-like surface based on spherical harmonics. Non-negative scalar 
%   alpha is a regularisation parameter and real scalar s is the parameter 
%   of the Sobolev seminorm of H^{s}(S, R).
%
%   X is an n-by-3 matrix, N are the degrees of spherical harmonics and
%   must be non-negative consecutive integers.
%
%   c is a vector of coefficients of the spherical harmonics Y and has
%   length dim. Y itself is an n-by-dim matrix of scalar spherical
%   harmonics evaluated at X projected to the unit sphere.
assert(size(X, 2) == 3);
assert(isscalar(alpha));
assert(alpha >= 0);
assert(isscalar(s));
assert(isvector(N));
assert(all(N >= 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Project data to unit sphere.
len = sqrt(sum(X.^2, 2));
X = bsxfun(@rdivide, X, len);

% Compute spherical harmonics at projected points.
Y = spharmn(N, X);
dim = size(Y, 2);

% Create symmetric system matrix.
A = zeros(dim, dim);
for p=1:dim
    for q=1:p
        A(p, q) = Y(:, p)'*Y(:, q);
        A(q, p) = A(p, q);
    end
end

% Create diagonal matrix.
D = alpha * diag(spharmeigs(N) .^ s);

% Compute vector b.
b = sum(bsxfun(@times, Y, len), 1)';

% Solve linear system.
c = gmres(A + D, b, [], 1e-6, min(1000, size(A, 1)));

end