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
function [Uf, Vf, u, v, L] = ofdbsolve(dim1, dim2, U, V, W, d1, d2, Y1, Y2, b, alpha, beta, s1, s2)
%OFDBSOLVE Solves the linear system where U and V are represented in
%different bases.
%
%   [Uf, Vf, u, v, L] = OFDBSOLVE(dim1, dim2, U, V, W, d1, d2, Y1, Y2, b, alpha, beta, s1, s2)
%   takes precomputed functions and solves the actual linear system for 
%   optical flow decomposition.
%
%   Uf and Vf are vector fields defined on the faces of the triangulation 
%   and are of size [n, 3], where n = size(Y1, 1) = size(Y2, 1) is the 
%   number of faces.
%
%   u and v are vectors of coefficients of vector spherical harmonics Y1
%   and Y2.
%
%   L is a struct containing information about the linear system solve.

assert(alpha >= 0);
assert(beta >= 0);
assert(dim1 > 0);
assert(dim2 > 0);
assert(isvector(d1));
assert(isvector(d2));
assert(length(d1) == dim1);
assert(length(d2) == dim2);
assert(size(U, 1) == dim1);
assert(size(U, 2) == dim1);
assert(size(V, 1) == dim2);
assert(size(V, 2) == dim2);
assert(size(W, 1) == dim1);
assert(size(W, 2) == dim2);
assert(size(Y1, 2) == dim1);
assert(size(Y2, 2) == dim2);
assert(size(Y1, 1) == size(Y2, 1));
assert(isscalar(s1));
assert(isscalar(s2));

% Get number of faces.
n = size(Y1, 1);

% Compute coefficients of norms.
ds = [alpha * (d1 .^ s1); beta * (d2 .^ s2)];

% Create function handle.
function v = fun(x)
    v1 = U * x(1:dim1) + W * x(dim1+1:dim1+dim2);
    v2 = W' * x(1:dim1) + V * x(dim1+1:dim1+dim2);
    v = [v1; v2] + ds .* x;
end

% Store norm of rhs.
L.rhs = norm(b, 2);

% Solve linear system.
ticId = tic;
[z, flag, relres, iter, resvec] = gmres(@fun, b, [], 1e-6, 30);

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
u = z(1:dim1);
v = z(dim1+1:end);

% Recover vector field.
Uf = zeros(n, 3);
parfor k=1:dim1
    Uf = Uf + u(k) * squeeze(Y1(:, k, :));
end
Vf = zeros(n, 3);
parfor k=1:dim2
    Vf = Vf + v(k) * squeeze(Y2(:, k, :));
end

end