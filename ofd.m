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
function [U, V] = ofd(N, f1, f2, h, alpha, beta)
%OFD Returns an optical flow decomposition.
%
%   [U, V] = OFD(N, f1, f2, h, alpha, beta) takes equally sized images f1
%   and f2, a spacing parameter h, and regularisation parameters alpha and
%   beta, returns an optical flow decomposition U, V.
%
%   U and V are of size [m, n, 3], where [m, n] = size(f1).

% Assumes static parametrisation.

[m, n] = size(f1);
assert(size(f2, 1) == m);
assert(size(f2, 2) == n);

% Compute time derivative of f.
dfdt = (f2 - f1) / h;

% Reproduce grid.
%[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));
t = linspace(-1, 1, n+2);
[phi, t] = ndgrid(linspace(-pi, pi, m), t(2:end-1));

% Vector containing eigenvalues.
D = zeros(N^2 + 2*N, 1);

Y1 = zeros(m, n, N^2 + 2*N, 3);
Y2 = zeros(m, n, N^2 + 2*N, 3);

% Create spherical harmonics.
c = 1;
for k=1:N
    [Yi, Yj] = vspharm(k, phi(:), t(:), m, n);
    Y1(:, :, c:(c+2*k), :) = Yi;
    Y2(:, :, c:(c+2*k), :) = Yj;
    % Save eigenvalues.
    D(c:(c+2*k)) = repmat(k*(k+1), (c+2*k)-c+1, 1);
    c = c + 2*k + 1;
end
Y = cat(3, Y1, Y2);
D = repmat(D, 2, 1);

%% Compute surface gradient of f1.
[dt, dphi] = gradient(f1, h);

% Compute orthonormal basis.
ephi(:, :, 1) = -sin(phi);
ephi(:, :, 2) = cos(phi);
ephi(:, :, 3) = 0;

et(:, :, 1) = -t .* cos(phi);
et(:, :, 2) = -t .* sin(phi);
et(:, :, 3) = sqrt(1 - t.^2);

gradf1(:, :, 1) = (ephi(:, :, 1) ./ sqrt(1 - t.^2)) .* dphi + sqrt(1 - t.^2) .* et(:, :, 1) .* dt;
gradf1(:, :, 2) = (ephi(:, :, 2) ./ sqrt(1 - t.^2)) .* dphi + sqrt(1 - t.^2) .* et(:, :, 2) .* dt;
gradf1(:, :, 3) = (ephi(:, :, 3) ./ sqrt(1 - t.^2)) .* dphi + sqrt(1 - t.^2) .* et(:, :, 3) .* dt;

%% Compute inner products grad f \cdot Y.
Z = zeros(m, n, 2*(N^2 + 2*N), 1);
for k=1:2*(N^2 + 2*N)
    Z(:, :, k) = gradf1(:, :, 1) .* Y(:, :, k, 1) + gradf1(:, :, 2) .* Y(:, :, k, 2) + gradf1(:, :, 3) .* Y(:, :, k, 3);
end

% Create matrix A tilde.
At = zeros(2*(N^2 + 2*N), 2*(N^2 + 2*N));
for p=1:2*(N^2 + 2*N)
    for q=1:p
        P = Z(:, :, p) .* Z(:, :, q);
        At(p, q) = sphericalIntegral(phi, t, P);
        At(q, p) = At(p, q);
    end
end

% Create vector b.
b = zeros(2*(N^2 + 2*N), 1);
for p=1:2*(N^2 + 2*N)
    P = dfdt .* Z(:, :, p);
    R = trapz(phi(:, 1), P, 1);
    b(p) = - trapz(t(1, :), R, 2);
end

% Create system matrix A.
A = [At + spdiags(alpha*D, 0, 2*(N^2 + 2*N), 2*(N^2 + 2*N)), At; At, At + spdiags(beta ./ D, 0, 2*(N^2 + 2*N), 2*(N^2 + 2*N))];
clear At;
clear Z;
clear D;

% Solve linear system.
z = cgs(A, repmat(b, 2, 1), 10e-6, 30);

% Recover coefficients.
u = z(1:2*(N^2 + 2*N));
v = z(2*(N^2 + 2*N)+1:end);

% Recover vector field.
U = zeros(m, n, 3);
V = zeros(m, n, 3);
for p=1:2*(N^2 + 2*N)
    U = U + u(p) * squeeze(Y(:, :, p, :));
    V = V + v(p) * squeeze(Y(:, :, p, :));
end

end