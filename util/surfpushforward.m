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
function Z = surfpushforward(N, c, F, V, Y)
%SURFPUSHFORWARD Computes the pushforward of tangent vector fields on the 
%unit sphere to a sphere-like surface.
%
%   Z = SURFPUSHFORWARD(N, c, F, V, Y) takes coefficients c for scalar 
%   sphercial harmonics of degrees N defining a sphere-like surface, a 
%   triangulation F, V of the unit sphere, and a piecewise constant vector 
%   field Y defined on the faces F and returns the pushforward Z of Y.
%
%   Note that the Y lives at the incenters of the faces F! Thus, also Z
%   lives on the incenters of the faces pushed to the surface.
%
%   Note that N must be a vector of non-negative consecutive integers. 
%   c is a vector of size dim, where dim is the number of scalar spherical 
%   harmonics of degrees specified by N. Y is an m-by-k-by-3 matrix 
%   defining vectors. F is an m-by-3 matrix and and V is an n-by-3 matrix 
%   of vertices on the unit sphere.
%
%   Z is a matrix of size m-by-k-by-3.

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N >= 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Compute and check dimension.
dim = N(end)^2 + 2*N(end) - N(1)^2 + 1;
assert(isvector(c));
assert(length(c) == dim);

% Check if dimensions comply.
assert(size(Y, 1) == size(F, 1));
assert(size(Y, 3) == 3);

% Compute triangle incenters and project to unit sphere.
TR = TriRep(F, V);
IC = TR.incenters;
len = sqrt(sum(IC .^ 2, 2));
IC = bsxfun(@rdivide, IC, len);

% Evaluate scalar spherical harmonics at incenters.
YIC = spharmn(N, IC);

% Evaluate scalar spherical harmonics at vertices.
YV = spharmn(N, V);

% Recover surface function at incenters.
rhoic = YIC * c;

% Recover surface function at vertices.
rhov = YV * c;

% Compute gradient of rho on triangulation.
g = grad(F, V, rhov);

% Compute pushforward for each vector field.
Z = zeros(size(Y));
for k=1:size(Y, 2)
    Yk = squeeze(Y(:, k, :));
    Z(:, k, :) = bsxfun(@times, Yk, rhoic) + bsxfun(@times, IC, sum(g .* Yk, 2));
end