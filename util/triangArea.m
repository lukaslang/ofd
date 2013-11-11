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
function A = triangArea(F, V)
%TRIANGAREA Computes the areas of the triangles for a given triangulation.
%
%   A = triangArea(F, V) takes a triangulation F, V and returns a vector A
%   containing the area of each triangle.
%
%   Note that size(A) = [size(F, 1), 1].

assert(size(F, 2) == 3);
assert(size(V, 2) == 3);

% Compute cross product of sides.
X = cross(V(F(:, 2), :) - V(F(:, 1), :), V(F(:, 3), :) - V(F(:, 1), :));

% Compute areas of triangles.
A = sqrt(sum(X .^2, 2)) ./ 2;

end