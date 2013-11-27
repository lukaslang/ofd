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
function [u, L] = ofsolve(dim, U, b, d, alpha, s)
%OFSOLVE Solves the linear system.
%
%   [u, L] = OFSOLVE(dim, U, b, d, alpha, s) takes precomputed 
%   functions and solves the linear system (U+diag(d))*u = b.
%
%   u is a vector of coefficients of the vector spherical harmonics basis 
%   and is of length dim.
%
%   L is a struct containing information about the linear system solve.

assert(alpha >= 0);
assert(dim > 0);
assert(isvector(d));
assert(length(d) == dim);
assert(all(size(U, 1) == [dim, dim]));
assert(isscalar(s));

% Compute coefficients of norm.
ds = alpha * (d .^ s);

% Create function handle.
function v = fun(x)
    v = U * x + ds .* x;
end

% Store norm of rhs.
L.rhs = norm(b, 2);

% Solve linear system.
ticId = tic;
[u, flag, relres, iter, resvec] = gmres(@fun, b, [], 1e-6, 30);

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

end