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
function v = sphericalIntegral(phi, t, f)
%SPHERICALINTEGRAL Numerically integrates a function on the unit sphere.
%
%   v = SPHERICALINTEGRAL(phi, t, f) takes a cylindrical parametrisation
%   (phi, t) each of size m-times-n and data f of same size and returns the
%   surface integral over f on the unit sphere.

I = trapz(phi(:, 1), f, 1);
v = trapz(t(1, :), I);

end