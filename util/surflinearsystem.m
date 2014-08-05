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
function [dim, A, D, b] = surflinearsystem(F, V, Ns, c, N, f1, f2, h, tol)
%SURFLINEARSYSTEM Computes the linear system used in optical flow on a
%surface.
%
%   [dim, A, D, b] = SURFLINEARSYSTEM(F, V, N, f1, f2, h, tol)
%   takes a triangulation F, V, degrees N, images f1, f2 on the vertices V 
%   of the triangulation, and a finite-difference parameter h for the time 
%   derivative of f1, f2 and returns the linear system needed for optical 
%   flow computation on the surface. The scalar tol is a tolerance 
%   parameter for numerical integration.
%
%   Note that degrees N must be a vector of consecutive positive integers!
%
%   The linear system returned is
%
%   (A + D) * x = b.
%
%   dim is the dimension of the linear system.
%
%   A amd D is a matrices of size dim-by-dim.
%
%   b is the right hand side and is of length dim.

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

assert(h > 0);
assert(isvector(f1));
assert(isvector(f2));
assert(~isempty(F));
assert(~isempty(V));
assert(~isempty(f1));
assert(~isempty(f2));

% Compute surface coordinates of vertices.
Vs = surfsynth(Ns, V, c);

% Compute approximate time derivative for each face.
dfdt = sum(f2(F) - f1(F), 2) ./ (3*h);

% Compute triangle areas to be used in integration.
a = triangArea(F, Vs);

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute surface gradient of first image.
gradf = grad(F, Vs, f1);

% Find indices where grad f is greater than epsilon. In areas where the
% length of the image gradient is almost zero inner products and thus
% surface integrals will be alomost zeros and thus can be excluded from
% numerical integration.
idx = sqrt(sum(gradf.^2, 2)) > tol;

% Constrain data.
Fc = F(idx, :);
gradfc = gradf(idx, :);
ac = a(idx);
dfdtc = dfdt(idx);

% Return all zero matrices if no data given.
if(numel(find(idx)) == 0)
    A = sparse([], [], [], dim, dim);
    D = A;
    b = sparse([], [], [], dim, 1);
    return;
end

% Compute inner products grad f with pushforward of vspharm.
Zc = surfvspharmdot(gradfc, Fc, V, N, Ns, c);

% Create matrix A.
A = matrixU(dim, Zc, Fc, V, ac);
    
% Compute orthonormal basis of the tangent space at incenters.
[B, ~] = surftangentialbasis(Ns, c, F, V);
[~, E1, E2] = orthonormalise(squeeze(B(:, 1, :)), squeeze(B(:, 2, :)));

% Compute covariant derivatives of vector spherical harmonics.
[D1E1, D1E2, D2E1, D2E2] = surfvspharmgrad(E1, E2, F, V, N, Ns, c);

% Create matrix D.
D = zeros(dim, dim);
for p=1:dim
    for q=1:p
        I = D1E1(:, p) .* D1E1(:, q) + D1E2(:, p) .* D1E2(:, q) + D2E1(:, p) .* D2E1(:, q) + D2E2(:, p) .* D2E2(:, q);
        D(p, q) = triangIntegral(F, Vs, I, a);
        D(q, p) = D(p, q);
    end
end

% Add diagonal.
d = zeros(dim, 1);
for p=1:dim
    d(p) = triangIntegral(F, Vs, D1E1(:, p).^2 + D1E2(:, p).^2 + D2E1(:, p).^2 + D2E2(:, p).^2, a);
end
D = diag(d) + D / 2;

% Create vector b.
b = zeros(dim, 1);
parfor k=1:dim
    b(k) = - triangIntegral(Fc, Vs, dfdtc .* Zc(:, k), ac);
end

end