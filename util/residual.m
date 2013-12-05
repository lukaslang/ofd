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
function [res, min] = residual(v, F, V, f1, f2, tol)
%RESIDUAL Computes the optical flow residual and the minimum residual on a 
%triangulation for a given velocity field.
%
%   [res, min] = residual(v, F, V, f1, f2, tol) takes a vector field v, a 
%   triangulation F, V, and images f1, f2 defined on the vertices V and 
%   returns the residual res of the optical flow equation and the minimum 
%   residual. tol is a tolerance parameter for numerical integration.
%
%   Note that f must be of size [n, 3], where n = size(F, 1) is the number
%   of faces. Note that res and min are scalars.

assert(all(size(v) == [size(F, 1), 3]));
assert(isvector(f1));
assert(isvector(f2));
assert(length(f1) == size(V, 1));
assert(length(f2) == size(V, 1));

% Compute approximate time derivative for each face.
dfdt = sum(f2(F) - f1(F), 2) ./ 3;

% Compute surface gradient of first image.
gradf = grad(F, V, f1);

% Find indices where grad f is greater than epsilon. In areas where the
% length of the image gradient is almost zero inner products and thus
% surface integrals will be alomost zeros and thus can be excluded from
% numerical integration.
idx = sqrt(sum(gradf.^2, 2)) > tol;

% Compute triangle areas.
a = triangArea(F, V);

% Compute (grad f \cdot v + dfdt).^2 for each face..
f = (dot(gradf(idx, :), v(idx, :), 2) + dfdt(idx)) .^2;

% Compute surface integral.
res = triangIntegral(F(idx, :), V, f, a(idx));

% Compute minimal residual.
min = triangIntegral(F(~idx), V, dfdt(~idx).^2, a(~idx));

end