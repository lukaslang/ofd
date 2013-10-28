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
function [x, y, z] = sphericalBand(m, n, r)
%SPHERICALBAND returns cartesian coordinates for a spherical band.
%
%   [x, y, z] = SPHERICALBAND(m, n, r) creates a m-times-n mesh in
%   spherical coordinates with different radii r, which is a vector, and 
%   returns cartesian coordinates.
%
%   Note that for a scalar r this function returns just the cartesian 
%   coordinates of a regular sphere of radius r!
%
%   Calling SPHERICALBAND(n, n, 1) returns basically the same as
%   SPHERE(n-1) except for the order of the vertices.

phi = linspace(-pi/2, pi/2, m);
theta = linspace(-pi, pi, n);

% Create spherical mesh.
[THETA, PHI, R] = meshgrid(phi, theta, r);

% Convert to cartesian coordinates.
[x, y, z] = sph2cart(PHI, THETA, R);

end