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
function [U, u] = ofsolve(dim, At, b, Y, d, alpha)
%OFSOLVE Solves the linear system.
%
%   [U, u] = OFSOLVE(dim, At, b, Y, d, alpha) takes precomputed functions 
%   and solves the actual linear system for optical flow.
%
%   U is defined on the faces of the triangulation and is of size [n, 3], 
%   where n = size(Y, 1) is the number of faces.
%
%   u is a vector of coefficients of vector spherical harmonics Y.

assert(alpha > 0);
assert(dim > 0);
assert(isvector(d));
assert(length(d) == dim);
assert(size(At, 1) == dim);
assert(size(At, 2) == dim);
assert(size(Y, 2) == dim);

% Get number of faces.
n = size(Y, 1);

% Create function handle.
function v = fun(x)
    v = At*x + alpha * d .* x;
end

% Solve linear system.
u = cgs(@fun, b, 10e-6, 30);

% Recover vector field.
U = zeros(n, 3);
parfor k=1:dim
    U = U + u(k) * squeeze(Y(:, k, :));
end

end