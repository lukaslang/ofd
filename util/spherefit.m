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
function [c, r] = spherefit(X, c, r)
%SPHEREFIT Fits a sphere to data.
%
%   [c, r] = SPHEREFIT(X) takes an n-by-3 matrix X and fits a sphere with 
%   center c and radius r. In addition, the function takes c and r as 
%   initial values.
%
%   Note that c is a 1-by-3 vector and r is a scalar.

assert(size(X, 2) == 3);
assert(size(c, 1) == 1);
assert(size(c, 2) == 3);
assert(isscalar(r));

% Define objective.
fun = @(x) sum(((X(:, 1) - x(1)).^2 + (X(:, 2) - x(2)).^2 + (X(:, 3) - x(3)).^2 - x(4).^2).^2);

% Find minimum.
x = fminsearch(fun, [c, r]);

% Recover solution.
c = x(1:3);
r = x(4);

end