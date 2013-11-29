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
function [dim1, dim2, U, V, W, d1, d2, b] = linearsystemdb(F, V, M, N, f1, f2, h, tol)
%LINEARSYSTEMDB Generates the computationally expensive parts of the 
%decomposition for different bases for later use.
%
%   [dim1, dim2, U, V, W, d1, d2, b] = LINEARSYSTEMDB(F, V, M, N, f1, f2, h, tol)
%   takes a triangulation F, V, degrees M and N, images f1, f2 on the 
%   vertices V of the triangulation, and a finite-difference parameter h 
%   for the time derivative of f1, f2 and returns the linear system needed 
%   for optical flow computation and the decomposition. The scalar tol is a
%   tolerance parameter for numerical integration.
%
%   Note that degrees M, N must be vectors of consecutive positive integers!
%
%   The linear system returned is
%   
%   [U+diag(d1),    V;
%    V',        W+diag(d2)] * [u; v] = b.
%
%   If M and N are the same intervals then all four matrices are the same 
%   and solving Ax=b can be implemented faster in an iterative manner. Use 
%   linearsystem instead!
%
%   dim1, dim2 are scalars so that dim1 + dim2 is the dimension of the linear
%   system.
%
%   d1, d2 are a vectors of length dim1 and dim2, respectively, and contain
%   the eigenvalues of the respective bases.
%
%   b is the right hand side and is of length dim1 + dim2.

% Check if both M and N are intervals of consecutive positive integers.
assert(isvector(M));
assert(all(M > 0));
assert(length(M) == M(end) - M(1) + 1);
assert(all((M == (M(1):M(end)))));
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
dfdt = sum(f2(F) - f1(F), 2) ./ 3;

% Compute triangle areas to be used in integration.
a = triangArea(F, V);

% Compute dimensions.
dim1 = 2*(M(end)^2 + 2*M(end) - M(1)^2 + 1);
dim2 = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

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
Z1 = vspharmdot(gradfc, Fc, V, M);
Z2 = vspharmdot(gradfc, Fc, V, N);

% Compute eigenvalues.
d1 = vspharmeigs(M);
d2 = vspharmeigs(N);

% Create matrices U and V.
U = matrixU(dim1, Z1, Fc, V, ac);
V = matrixU(dim2, Z2, Fc, V, ac);

% Create matrix W.
W = zeros(dim1, dim2);
for p=1:dim1
    for q=1:dim2
        W(p, q) = triangIntegral(Fc, V, Z1(:, p) .* Z2(:, q), ac);
    end
end

% Create vector b.
b = zeros(dim1 + dim2, 1);
parfor k=1:dim1
    b(k) = - triangIntegral(Fc, V, dfdtc .* Z1(:, k), ac);
end
parfor k=1:dim2
    b(dim1 + k) = - triangIntegral(Fc, V, dfdtc .* Z2(:, k), ac);
end

end