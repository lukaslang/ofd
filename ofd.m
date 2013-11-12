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

% Vector containing eigenvalues.
D = zeros(N^2 + 2*N, 1);

Y1 = zeros(n, N^2 + 2*N, 3);
Y2 = zeros(n, N^2 + 2*N, 3);

% Create spherical harmonics.
c = 1;
for k=1:N
    [Yi, Yj] = vspharm(k, F, V);
    Y1(:, c:(c+2*k), :) = Yi;
    Y2(:, c:(c+2*k), :) = Yj;
    % Save eigenvalues.
    D(c:(c+2*k)) = repmat(k*(k+1), (c+2*k)-c+1, 1);
    c = c + 2*k + 1;
end
Y = cat(2, Y1, Y2);
D = repmat(D, 2, 1);

% Compute surface gradient of first image.
gradf = grad(F, V, f1);

% Compute inner products grad f \cdot Y.
disp('Computing inner products grad f cdot Y.');
Z = zeros(n, 2*(N^2 + 2*N));
tic;
for k=1:2*(N^2 + 2*N)
    Z(:, k) = dot(gradf, squeeze(Y(:, k, :)), 2);
end
toc;

% Create matrix A tilde.
disp('Computing A tilde.');
At = zeros(2*(N^2 + 2*N), 2*(N^2 + 2*N));
tic;
for p=1:2*(N^2 + 2*N)
    for q=1:p
        At(p, q) = triangIntegral(F, V, Z(:, p) .* Z(:, q), a);
        At(q, p) = At(p, q);
    end
end
toc;
clear P;

% Create vector b.
disp('Computing vector b.');
b = zeros(2*(N^2 + 2*N), 1);
tic;
for p=1:2*(N^2 + 2*N)
    b(p) = - triangIntegral(F, V, dfdt .* Z(:, p), a);
end
toc;
clear Z;

% Create system matrix A.
disp('Solving linear system...');
tic;
A = [At + spdiags(alpha*D, 0, 2*(N^2 + 2*N), 2*(N^2 + 2*N)), At; At, At + spdiags(beta ./ D, 0, 2*(N^2 + 2*N), 2*(N^2 + 2*N))];
toc;
clear At;
clear D;

% Solve linear system.
z = cgs(A, repmat(b, 2, 1), 10e-6, 30);

% Recover coefficients.
u = z(1:2*(N^2 + 2*N));
v = z(2*(N^2 + 2*N)+1:end);

% Recover vector field.
disp('Recover vector field.');
U = squeeze(sum(repmat(u', [n, 1, 3]) .* Y, 2));
V = squeeze(sum(repmat(v', [n, 1, 3]) .* Y, 2));

end