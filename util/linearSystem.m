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
function [dim, U, d, b] = linearsystem(F, V, N, f1, f2, h, tol)
%LINEARSYSTEM Computes the linear system used in optical flow and optical
%flow decomposition.
%
%   [dim, U, d, b] = LINEARSYSTEM(F, V, N, f1, f2, h, tol)
%   takes a triangulation F, V, degrees N, images f1, f2 on the vertices V 
%   of the triangulation, and a finite-difference parameter h for the time 
%   derivative of f1, f2 and returns the linear system needed for optical 
%   flow computation and the decomposition. The scalar tol is a tolerance 
%   parameter for numerical integration.
%
%   Note that degrees N must be a vector of consecutive positive integers!
%
%   The linear system returned is
%
%   (U + diag(d)) * x = b.
%
%   dim is the dimension of the linear system.
%
%   d is a vector of length dim and contains the eigenvalues of the basis.
%
%   b is the right hand side and is of length dim.

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

assert(h > 0);
assert(isvector(f1));
assert(isvector(f2));
assert(~isempty(F));
assert(~isempty(V));
assert(~isempty(f1));
assert(~isempty(f2));

% Compute approximate time derivative for each face.
dfdt = sum(f2(F) - f1(F), 2) ./ (3*h);

% Compute triangle areas to be used in integration.
a = triangArea(F, V);

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute surface gradient of first image.
gradf = grad(F, V, f1);

% Find indices where grad f is greater than epsilon. In areas where the
% length of the image gradient is almost zero inner products and thus
% surface integrals will be alomost zeros and thus can be excluded from
% numerical integration.
idx = sqrt(sum(gradf.^2, 2)) > tol;

% Constrain data.
Fc = F(idx, :);
gradfc = gradf(idx, :);
ac = a(idx);
dfdtc = dfdt(idx);

% Compute inner products grad f \cdot Y.
Z = vspharmdot(gradfc, Fc, V, N);

% Compute eigenvalues.
d = vspharmeigs(N);

% Create matrix U.
U = matrixU(dim, Z, Fc, V, ac);

% Create vector b.
b = zeros(dim, 1);
parfor k=1:dim
    b(k) = - triangIntegral(Fc, V, dfdtc .* Z(:, k), ac);
end

end