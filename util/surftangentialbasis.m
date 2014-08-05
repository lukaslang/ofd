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
function [Z, IC] = surftangentialbasis(Ns, c, F, V)
%SURFTANGENTIALBASIS Computes a tangential basis on a sphere-like surface.
%
%   [Z, IC] = SURFTANGENTIALBASIS(N, c, F, V) takes coefficients c for 
%   scalar sphercial harmonics of degrees Ns defining a sphere-like surface,
%   a triangulation F, V of the unit sphere and returns a tangential basis
%   Z at the incenters IC of the faces on the surface.
%
%   NOTE: the basis Z is computed at the incenters and therefore might NOT
%   be tangential to the faces on the surface!
%
%   Note that Ns must be a vector of non-negative consecutive integers. 
%   c is a vector of size dim, where dim is the number of scalar spherical 
%   harmonics of degrees specified by Ns. F is an m-by-3 matrix and and V is
%   an n-by-3 matrix of vertices on the unit sphere.
%
%   Z is a matrix of size m-by-2-by-3. The tangential basis is given by
%   Z(:, i, :) for i={1, 2}. IC is a matrix of size m-by-3.

% Check if Ns is an interval of consecutive positive integers.
assert(isvector(Ns));
assert(all(Ns >= 0));
assert(length(Ns) == Ns(end) - Ns(1) + 1);
assert(all((Ns == (Ns(1):Ns(end)))));

% Compute and check dimension.
dim = Ns(end)^2 + 2*Ns(end) - Ns(1)^2 + 1;
assert(isvector(c));
assert(length(c) == dim);

% Compute triangle incenters and project to unit sphere.
TR = TriRep(F, V);
IC = TR.incenters;
len = sqrt(sum(IC .^ 2, 2));
IC = bsxfun(@rdivide, IC, len);

% Convert to cylindrical coordinates.
[phi, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
t = sin(pi - el);

% Compute orthonormal basis on the unit sphere.
B(:, 1, 1) = -sin(phi);
B(:, 1, 2) = cos(phi);
B(:, 1, 3) = 0;
B(:, 2, 1) = -t .* cos(phi);
B(:, 2, 2) = -t .* sin(phi);
B(:, 2, 3) = sqrt(1 - t.^2);

% Compute pushforward of spherical basis and incenters.
[Z, IC] = surfpushforward(Ns, c, F, V, B);