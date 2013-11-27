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
function Fn = facenormals(F, V)
%FACENORMALS Computes the face unit normals for a given triangulation.
%
%   Fn = FACENORMALS(F, V) takes a triangulation F, V and returns a matrix
%   Fn of face unit normals. Triangles are assumed to be oriented clockwise.
%
%   Note that size(Fn) = [size(F, 1), 3].

assert(size(F, 2) == 3);
assert(size(V, 2) == 3);

% Compute cross product of sides.
Fn = cross(V(F(:, 3), :) - V(F(:, 1), :), V(F(:, 2), :) - V(F(:, 1), :));

% Normalise.
Fn = bsxfun(@rdivide, Fn, sqrt(sum(Fn .^ 2, 2)));

end