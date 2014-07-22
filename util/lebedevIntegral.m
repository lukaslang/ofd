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
function v = lebedevIntegral(f, deg)
%LEBEDEVINTEGRAL Applies the Lebedev cubature to numerically integrate a
%a function on the unit sphere.
%
%   v = LEBEDEVINTEGRAL(f, deg) takes the handle of a function f defined 
%   on the unit sphere and returns the surface integral int_S f(x, y, z) dS.
%   f is assumed to be a function which can be expanded in scalar spherical
%   harmonics of degree up to deg. The number of evaluation points of the 
%   Lebedev cubature is implied by the given degree.
%
%   The function f is required to have the signature fv = f(x, y, z), where
%   all arguments are vectors of length n, (x, y, z) are arbitrary points 
%   on the unit sphere. Note that fv must also be a vector of length n.
%
%   Example: f = @(x, y, z) x.^2 + y.^2 + z.^2 is constant 1.
assert(isa(f, 'function_handle'));
assert(isscalar(deg));
assert(deg >= 0);

% Implemented degrees.
degs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41,...
    47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131];
% Associated number of points.
pts = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074,...
    3470, 3890, 4334, 4802, 5294, 5810];

% Check if degree is supported by Lebedev cubature.
[~, idx] = ismember(deg, degs, 'R2012a');
assert(isscalar(idx));
assert(idx > 0);

% Get pre-computed Lebedev cubature points.
leb = getLebedevSphere(pts(idx));

% Evaluate function f at cubature points.
fv = f(leb.x, leb.y, leb.z);

% Compute integral approximation.
v = fv' * leb.w;

end