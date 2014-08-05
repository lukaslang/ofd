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
function [Z, ICs] = surfpushforward(Ns, c, F, V, Y)
%SURFPUSHFORWARD Computes the pushforward of tangent vector fields on the 
%unit sphere to a sphere-like surface.
%
%   [Z, ICs] = SURFPUSHFORWARD(Ns, c, F, V, Y) takes coefficients c for scalar 
%   sphercial harmonics of degrees Ns defining a sphere-like surface, a 
%   triangulation F, V of the unit sphere, and a piecewise constant vector 
%   field Y defined on the faces F and returns the pushforward Z of Y. In
%   addition, the pushforward ICs of the incenters to the surface are 
%   returned.
%
%   Note that the Y lives at the incenters of the faces F! Thus, also Z
%   lives on the incenters of the faces pushed to the surface.
%
%   Note that Ns must be a vector of non-negative consecutive integers. 
%   c is a vector of size dim, where dim is the number of scalar spherical 
%   harmonics of degrees specified by Ns. Y is an m-by-k-by-3 matrix 
%   defining vectors. F is an m-by-3 matrix and and V is an n-by-3 matrix 
%   of vertices on the unit sphere.
%
%   Z is a matrix of size m-by-k-by-3. ICs is an m-by-3 matrix.

% Check if N is an interval of consecutive positive integers.
assert(isvector(Ns));
assert(all(Ns >= 0));
assert(length(Ns) == Ns(end) - Ns(1) + 1);
assert(all((Ns == (Ns(1):Ns(end)))));

% Compute and check dimension.
dim = Ns(end)^2 + 2*Ns(end) - Ns(1)^2 + 1;
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

% Compute surface coordinates of incenters.
[ICs, rhoic] = surfsynth(Ns, IC, c);

% Compute surface coordinates of vertices.
[~, rhov] = surfsynth(Ns, V, c);

% Compute gradient of rho on triangulation.
g = grad(F, V, rhov);

% Compute pushforward.
Z = pushforward(Y, IC, rhoic, g);