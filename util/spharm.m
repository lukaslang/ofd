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
function Y = spharm(N, X)
%SPHARM Generates fully normalised spherical harmonics for given cartesian
%coordinates.
%
%   Y = SPHARM(N, X) takes a n-by-3 matrix X of cartesian coordinates and 
%   returns fully normalised spherical harmonics Y_Nj of degree N >= 0 and 
%   j=-N,...,N. Note that every point X(k, :) must be on the unit sphere!
%
%   Note that size(Y) = [n, 2*N + 1].

assert(N >= 0);
assert(size(X, 2) == 3);

% Convert to cylindrical coordinates.
[phi, el, ~] = cart2sph(X(:, 1), X(:, 2), X(:, 3));
t = sin(pi - el);

% Compute fully normalised spherical harmonics.
Y = spharmp(N, phi, t);

end