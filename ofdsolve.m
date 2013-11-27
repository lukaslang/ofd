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
function [u, v, L] = ofdsolve(dim, U, b, d, alpha, beta, s1, s2)
%OFDSOLVE Solves the linear system.
%
%   [u, v, L] = OFDSOLVE(dim, U, b, d, alpha, beta, s1, s2) takes precomputed 
%   functions and solves the linear system
%
%   [U + diag(d), U; U', U + diag(d)]*[u; v] = [b; b].
%
%   u and v are vectors of coefficients of vector spherical harmonics basis.
%
%   L is a struct containing information about the linear system solve.

assert(alpha >= 0);
assert(beta >= 0);
assert(dim > 0);
assert(isvector(d));
assert(length(d) == dim);
assert(all(size(U) == [dim, dim]));
assert(isscalar(s1));
assert(isscalar(s2));

% Compute coefficients of norms.
ds = [alpha * (d .^ s1); beta * (d .^ s2)];

% Create function handle.
function v = fun(x)
    tv = U * (x(1:dim) + x(dim+1:end));
    v = repmat(tv, 2, 1) + ds .* x;
end

% Store norm of rhs.
rhs = repmat(b, 2, 1);
L.rhs = norm(rhs, 2);

% Solve linear system.
ticId = tic;
[z, flag, relres, iter, resvec] = gmres(@fun, rhs, [], 1e-6, 30);

% Store solver information.
L.time = toc(ticId);
L.flag = flag;
L.relres = relres;
L.iter = iter;
L.resvec = resvec;
L.tol = 1e-6;
L.maxit = 30;
L.solver = 'gmres';
L.restart = 0;

% Recover coefficients.
u = z(1:dim);
v = z(dim+1:end);

end