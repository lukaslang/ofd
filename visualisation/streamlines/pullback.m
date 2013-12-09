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
function w = pullback(P, v)
%PULLBACK Pulls back a vector field on a hemisphere to the coordinate
%domain by assuming the surface to be parametrised by a monge patch.
%
%   w = pullback(P, v) takes a vector field v at points P on a hemisphere
%   and pulls back v to the plane so that v = J*w.

% Assert that no point lies outside the unit circle.
assert(all(P(:, 1).^2 + P(:, 2).^2 <= 1));

% Compute third coordinate.
z = sqrt(1 - P(:, 1).^2 - P(:, 2).^2);

% Compute pullback of vector field.
PB(:, 1, 1) = 1 + P(:, 2).^2;
PB(:, 1, 2) = - P(:, 1) .* P(:, 2);
PB(:, 1, 3) = - P(:, 1) .* z;
PB(:, 2, 1) = - P(:, 1) .* P(:, 2);
PB(:, 2, 2) = 1 + P(:, 1).^2;
PB(:, 2, 3) = - P(:, 2) .* z;

% Get vector field w = g^{-1}*J'*v in the plane.
w(:, 1) = PB(:, 1, 1) .* v(:, 1) + PB(:, 1, 2) .* v(:, 2) + PB(:, 1, 3) .* v(:, 3);
w(:, 2) = PB(:, 2, 1) .* v(:, 1) + PB(:, 2, 2) .* v(:, 2) + PB(:, 2, 3) .* v(:, 3);

end