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
function Y = spharmp(N, phi, t)
%SPHARMP Generates fully normalised spherical harmonics on cylindrical
%parametrisation.
%
%   Y = SPHARMP(N, phi, t) takes matrices phi and t from a cylindrical 
%   parametrisation (e.g. created by ndgrid) and returns spherical 
%   harmonics Y_Nj of degree N >= 0 and j=-N,...,N. Note that phi must be 
%   in [0, 2pi[ and t in [-1, 1].
%
%   Note that size(Y) = [m*n, 2*N + 1].

[m, n] = size(phi);

assert(N >= 0);
assert(size(t, 1) == m);
assert(size(t, 2) == n);

% Compute legendre polynomials.
P = legendre(N, t(:), 'norm')';

% Fully normalise.
P = P .* [repmat(sqrt(2), m*n, 1), repmat(2, m*n, N)] ./ sqrt(4*pi);

% Reorder so that order is -n,...,n.
y = repmat(-N:N, m*n, 1) .* repmat(phi(:), 1, 2*N + 1);
y = [cos(y(:, 1:N + 1)), sin(y(:, N + 2:2*N + 1))];
Y = y .* [flip(P, 2), P(:, 2:end)];

end