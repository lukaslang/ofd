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
function [U, V, u, v, L] = ofdsolve(dim, At, b, Y, d, alpha, beta)
%OFDSOLVE Solves the linear system.
%
%   [U, V, u, v, L] = OFDSOLVE(dim, At, b, Y, d, alpha, beta) takes 
%   precomputed functions and solves the actual linear system for optical 
%   flow.
%
%   U and V are defined on the faces of the triangulation and is of size 
%   [n, 3], where n = size(Y, 1) is the number of faces.
%
%   u and v are vectors of coefficients of vector spherical harmonics Y.
%
%   L is a struct containing information about the linear system solve.

assert(alpha > 0);
assert(beta > 0);
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
    tv = At * (x(1:dim) + x(dim+1:end));
    v = repmat(tv, 2, 1) + [alpha .* d; beta ./ d] .* x;
end

% Store norm of rhs.
rhs = repmat(b, 2, 1);
L.rhs = norm(rhs, 2);

% Solve linear system.
ticId = tic;
[z, flag, relres, iter, resvec] = cgs(@fun, rhs, 1e-6, 30);

% Store solver information.
L.time = toc(ticId);
L.flag = flag;
L.relres = relres;
L.iter = iter;
L.resvec = resvec;
L.tol = 1e-6;
L.maxit = 30;

% Recover coefficients.
u = z(1:dim);
v = z(dim+1:end);

% Recover vector field.
U = zeros(n, 3);
V = zeros(n, 3);
parfor k=1:dim
    U = U + u(k) * squeeze(Y(:, k, :));
    V = V + v(k) * squeeze(Y(:, k, :));
end

end