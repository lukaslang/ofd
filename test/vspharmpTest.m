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
function test_suite = vspharmpTest
    initTestSuite;
end

function resultTest

% Create parametrisation and coordinates.
m = 10;
n = 20;

% Create parametrisation and polar coordinates.
[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));

% Create spherical harmonics.
N = 5;
[Y1, Y2] = vspharmp(N, phi, t);
assertFalse(isempty(Y1));
assertEqual(size(Y1), [m, n, 2*N + 1, 3]);
assertFalse(isempty(Y2));
assertEqual(size(Y2), [m, n, 2*N + 1, 3]);

end

function orthogonalityTest

% Create parametrisation and coordinates.
m = 10;
n = 20;

% Create parametrisation and polar coordinates.
t = linspace(-1, 1, n+2);
[phi, t] = ndgrid(linspace(-pi, pi, m), t(2:end-1));

% Create spherical harmonics.
N = 5;
[Y1, Y2] = vspharmp(N, phi, t);

% Compute R3 inner product.
ip = Y1(:, :, :, 1) .* Y2(:, :, :, 1) + Y1(:, :, :, 2) .* Y2(:, :, :, 2) + Y1(:, :, :, 3) .* Y2(:, :, :, 3);
assertAlmostEqual(ip, zeros(m, n, 2*N + 1));

end

function visualiseTest

% Create parametrisation and coordinates.
m = 50;
n = 50;

% Create parametrisation and polar coordinates.
t = linspace(-1, 1, n+2);
[phi, t] = ndgrid(linspace(-pi, pi, m), t(2:end-1));

% Create spherical harmonics.
N = 3;
[Y1, Y2] = vspharmp(N, phi, t);
assertFalse(isempty(Y1));
assertEqual(size(Y1), [m, n, 2*N + 1, 3]);
assertFalse(isempty(Y2));
assertEqual(size(Y2), [m, n, 2*N + 1, 3]);

% Convert to cartesian coordinates.
x = sqrt(1 - t .^2) .* cos(phi);
y = sqrt(1 - t .^2) .* sin(phi);
z = t;

% Create spherical harmonics for visualisation.
Ynj = spharmp(N, phi(:), t(:));

figure;
for k=1:2*N+1
    f = Ynj(:, k);
    f = reshape(f, m, n);

    % Visualise first kind.
    subplot(2*N+1, 2, 2*(k-1)+1);
    hold on;
    surf(x, y, z, f, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    quiver3(x, y, z, Y1(:, :, k, 1), Y1(:, :, k, 2), Y1(:, :, k, 3), 0);
    daspect([1, 1, 1]);
    axis([-1, 1, -1, 1, -1, 1]);
    view(3);

    % Visualise second kind.
    subplot(2*N+1, 2, 2*k);
    hold on;
    surf(x, y, z, f, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    quiver3(x, y, z, Y2(:, :, k, 1), Y2(:, :, k, 2), Y2(:, :, k, 3), 0);
    daspect([1, 1, 1]);
    axis([-1, 1, -1, 1, -1, 1]);
    view(3);
    
end
end