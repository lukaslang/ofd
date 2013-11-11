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
function v = triangIntegral(F, V, f)
%TRIANGINTEGRAL Takes a triangulation and numerically integrates a 
%function given on the faces.
%
%   v = TRIANGINTEGRAL(F, V, f) takes a triangulation F, V and a function f
%   defined on the faces F and returns the surface integral over f.

assert(isvector(f));
assert(size(f, 1) == size(F, 1));
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);

% Compute areas of triangles.
A = triangArea(F, V);

% Compute integral approximation.
v = A'*f;

end