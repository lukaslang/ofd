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
function Y = spharm(n, phi, t)
%SPHARM Generates fully normalised spherical harmonics.
%
%   Y = spharm(n, phi, t) takes vectors [phi, t] and returns spherical 
%   harmonics Y_nj of degree n >= 0 and j=-n,...,n. phi must be in [0, 2pi[ and
%   t in [-1, 1].
%
%   Note that size(Y) = [length(phi), 2n + 1].

m = length(phi);

assert(n >= 0);
assert(size(phi, 1) == m);
assert(size(phi, 2) == 1);
assert(size(t, 1) == m);
assert(size(t, 2) == 1);

% Compute legendre polynomials.
P = legendre(n, t, 'norm')';

% Fully normalise.
P = P .* [repmat(sqrt(2), m, 1), repmat(2, m, n)] ./ sqrt(4*pi);

% Reorder so that order is -n,...,n.
y = repmat(-n:n, m, 1) .* repmat(phi, 1, 2*n + 1);
y = [cos(y(:, 1:n+1)), sin(y(:, n+2:2*n+1))];
Y = y .* [flipdim(P, 2), P(:, 2:end)];

end