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
function [U, V] = ofd(N, F, V, f1, f2, h, alpha, beta)
%OFD Returns an optical flow decomposition on a triangulation.
%
%   [U, V] = OFD(N, F, V, f1, f2, h, alpha, beta) takes a triangulation 
%   F, V and images f1, f2 on the vertices of the triangulation and returns
%   an optical flow decomposition U, V. Scalar h is a spacing parameter and
%   alpha, beta are regularisation parameters.
%
%   U and V are defined on the faces F and are of size [n, 3], where 
%   n = size(F, 1) is the number of faces.

m = size(V, 1);
n = size(F, 1);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(isvector(f1));
assert(isvector(f2));
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);
assert(alpha > 0);
assert(beta > 0);
assert(h > 0);

% Compute approximate time derivative for each face.
dfdt = sum(f2(F) - f1(F), 2) ./ 3;

% Compute triangle areas to be used in integration.
a = triangArea(F, V);

% Create vector spherical harmonics up to degree N.
[Y, d] = vspharmn(N, F, V);

% Compute dimension.
dim = 2*(N^2 + 2*N);

% Compute surface gradient of first image.
gradf = grad(F, V, f1);

% Find indices where grad f is greater than epsilon. In areas where the
% length of the image gradient is almost zero inner products and thus
% surface integrals will be alomost zeros and thus can be excluded from
% numerical integration.
tol = 1e-6;
idx = sqrt(sum(gradf.^2, 2)) > tol;

% Constrain data.
Fc = F(idx, :);
nc = size(Fc, 1);
gradfc = gradf(idx, :);
Yc = Y(idx, :, :);
ac = a(idx);
dfdtc = dfdt(idx);

% Compute inner products grad f \cdot Y.
disp('Computing inner products grad f cdot Y.');
Z = zeros(nc, dim);
tic;
parfor k=1:dim
    Z(:, k) = dot(gradfc, squeeze(Yc(:, k, :)), 2);
end
toc;

% Create matrix A tilde.
disp('Computing A tilde.');
tic;
At = matrixAt(dim, Z, Fc, V, ac);
toc;

% Create vector b.
disp('Computing vector b.');
b = zeros(dim, 1);
tic;
parfor k=1:dim
    b(k) = - triangIntegral(Fc, V, dfdtc .* Z(:, k), ac);
end
toc;
clear Z;

% Create function handle.
function v = fun(x)
    tv = At * (x(1:dim) + x(dim+1:end));
    v = repmat(tv, 2, 1) + [alpha .* d; beta ./ d] .* x;
end

% Solve linear system.
disp('Solve linear system...');
tic;
z = cgs(@fun, repmat(b, 2, 1), 10e-6, 30);
toc;
clear b;

% Recover coefficients.
u = z(1:dim);
v = z(dim+1:end);

% Recover vector field.
disp('Recover vector field.');
U = zeros(n, 3);
V = zeros(n, 3);
parfor k=1:dim
    U = U + u(k) * squeeze(Y(:, k, :));
    V = V + v(k) * squeeze(Y(:, k, :));
end

end