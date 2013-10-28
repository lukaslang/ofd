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
function f = maximumIntensity(c, m, n, rs, X, Y, Z, u)
%MAXIMUMINTENSITY Returns the maximum intensity along radial lines.
%
%   f = MAXIMUMINTENSITY(c, m, n, rs, X, Y, Z, u) takes the center c of a spherical
%   band with radii rs, creates an m-times-n mesh and extracts the radial
%   maximum intensity from the data u which is defined on X, Y, and Z.

%   Note that size(f) is [m, n].

% Create spherical band.
[x, y, z] = sphericalBand(m, n, rs);

% Get data from cube.
f = dataFromCube(c(1)+x(:), c(2)+y(:), c(3)+z(:), X, Y, Z, u);

% Reshape and compute maximum intensity along third dimension.
f = max(reshape(f, m, n, length(rs)), [], 3);

end