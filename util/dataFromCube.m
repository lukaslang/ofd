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
function f = dataFromCube(x, y, z, X, Y, Z, u)
%DATAFROMCUBE Extracts data from a cube.
%
%   f = DATAFROMCUBE(x, y, z, X, Y, Z, u) takes points [x, y, z] and 
%   returns the interpolation of the function u, which is defined on 
%   [X, Y, Z] if exists and zero otherwise. X, Y, and Z usually arise from
%   calling ndgrid.
%
%   Note taht f is a vector of same length as x, y, and z.

assert(isvector(x));
assert(isequal(size(x), size(y), size(z)));

% Interpolate.
F = griddedInterpolant(X, Y, Z, u);
f = F(x, y, z);
f(isnan(f)) = 0;

end