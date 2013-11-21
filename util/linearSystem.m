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
function [dim1, dim2, U, V, W, d1, d2, Y1, Y2, b] = linearSystem(F, V, M, N, f1, f2, h, tol)
%LINEARSYSTEM Generates the computationally expensive parts of the 
%decomposition for different bases for later use.
%
%   [dim1, dim2, U, V, W, d1, d2, Y1, Y2, b] = LINEARSYSTEM(F, V, M, N, f1, f2, h, tol)
%   takes a triangulation F, V, degrees M and N, images f1, f2 on the 
%   vertices V of the triangulation, and a finite-difference parameter h 
%   for the time derivative of f1, f2 and returns the linear system needed 
%   for optical flow computation and the decomposition. The scalar tol is a
%   tolerance parameter for numerical integration.
%
%   Note that degrees must be either scalars M, N > 0 or vectors of 
%   consecutive positive integers!
%
%   The linear system returned is [U, V; V', W]*x = b. If M==N then all
%   four matrices are the same and are generated only once. Solving Ax=b 
%   can be implemented faster in an iterative manner. See
%   computeDataFunctions!
%
%   dim1, dim2 are scalars so that dim1 + dim2 is the dimension of the linear
%   system.
%
%   d1, d2 are a vectors of length dim1 and dim2, respectively, and contain
%   the eigenvalues of the bases Y1 and Y2.
%
%   b is the right hand side and is of length dim1 + dim2.

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

% Create vector spherical harmonics.
[Y1, d1] = vspharmn(M, F, V);
[Y2, d2] = vspharmn(N, F, V);

% Compute dimension.
dim1 = length(d1);
dim2 = length(d2);

% Compute surface gradient of first image.
gradf = grad(F, V, f1);

% Find indices where grad f is greater than epsilon. In areas where the
% length of the image gradient is almost zero inner products and thus
% surface integrals will be alomost zeros and thus can be excluded from
% numerical integration.
idx = sqrt(sum(gradf.^2, 2)) > tol;

% Constrain data.
Fc = F(idx, :);
nc = size(Fc, 1);
gradfc = gradf(idx, :);
Y1c = Y1(idx, :, :);
Y2c = Y2(idx, :, :);
ac = a(idx);
dfdtc = dfdt(idx);

% Compute inner products grad f \cdot Y.
Z1 = zeros(nc, dim1);
parfor k=1:dim1
    Z1(:, k) = dot(gradfc, squeeze(Y1c(:, k, :)), 2);
end
Z2 = zeros(nc, dim2);
parfor k=1:dim2
    Z2(:, k) = dot(gradfc, squeeze(Y2c(:, k, :)), 2);
end

% Create matrices U and V.
U = matrixAt(dim1, Z1, Fc, V, ac);
V = matrixAt(dim2, Z2, Fc, V, ac);

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