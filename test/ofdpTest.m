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
function test_suite = ofdpTest
    initTestSuite;
end

function resultTest

m = 30;
n = 20;

% Create two images.
f1 = zeros(m, n);
f2 = zeros(m, n);

N = 5;
h = 1;
alpha = 1;
beta = 1;

[u, v] = ofdp(N, f1, f2, h, alpha, beta);

assertFalse(isempty(u));
assertFalse(isempty(v));
assertEqual(size(u), [m, n, 3]);
assertEqual(size(v), [m, n, 3]);
assertEqual(u, zeros(m, n, 3));
assertEqual(v, zeros(m, n, 3));

end

function visualiseTest

% Create parametrisation and coordinates.
m = 100;
n = 100;

% Create parametrisation (without poles) and polar coordinates.
t = linspace(-1, 1, n+2);
[phi, t] = ndgrid(linspace(-pi, pi, m), t(2:end-1));

% Create two images.
Ynj = spharmp(5, phi, t);
f1 = 255 * reshape(Ynj(:, 3), m, n);
f2 = circshift(f1, [m/10, 0]);

N = 3;
h = 1;
alpha = 10;
beta = 1000;

[u, v] = ofdp(N, f1, f2, h, alpha, beta);

assertFalse(isempty(u));
assertFalse(isempty(v));
assertEqual(size(u), [m, n, 3]);
assertEqual(size(v), [m, n, 3]);

% Convert to cartesian coordinates.
x = sqrt(1 - t .^2) .* cos(phi);
y = sqrt(1 - t .^2) .* sin(phi);
z = t;

figure;
subplot(1, 2, 1);
surf(x, y, z, f1, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);
subplot(1, 2, 2);
surf(x, y, z, f2, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

figure;
subplot(1, 3, 1);
hold on;
surf(x, y, z, f1, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
quiver3(x, y, z, u(:, :, 1), u(:, :, 2), u(:, :, 3), 0, 'g');
daspect([1, 1, 1]);
view(3);

subplot(1, 3, 2);
hold on;
surf(x, y, z, f1, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
quiver3(x, y, z, v(:, :, 1), v(:, :, 2), v(:, :, 3), 0, 'y');
daspect([1, 1, 1]);
view(3);

subplot(1, 3, 3);
hold on;
surf(x, y, z, f1, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
quiver3(x, y, z, u(:, :, 1)+v(:, :, 1), u(:, :, 2)+v(:, :, 2), u(:, :, 3)+v(:, :, 3), 0, 'm');
daspect([1, 1, 1]);
view(3);


end