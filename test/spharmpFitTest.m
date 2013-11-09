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
function test_suite = spharmFitTest
    initTestSuite;
end

function resultTest

% Create parametrisation and coordinates.
m = 100;
n = 100;
r = 1;

% Create parametrisation and polar coordinates.
[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));

% Convert to cartesian coordinates.
x = r .* sqrt(1 - t .^2) .* cos(phi);
y = r .* sqrt(1 - t .^2) .* sin(phi);
z = r .* t;

% Create data.
[X, Y, Z] = ndgrid(-1:0.1:1, -1:0.1:1, -1:0.1:1);
u = X .* Y .* Z;
f = dataFromCube(x(:), y(:), z(:), X, Y, Z, u);
f = reshape(f, m, n);

% Visualise data.
figure;
subplot(1, 3, 1);
surf(x, y, z, f, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
colorbar;
daspect([1, 1, 1]);
view(3);

% Fit spherical harmonics.
N = 3;
[a, Ynj] = spharmpFit(N, phi, t, f);
assertFalse(isempty(a));
assertFalse(isempty(Ynj));
assertEqual(length(a), (N+1)^2);
assertEqual(size(Ynj), [m*n, (N+1)^2]);

% Recover fitted solution.
g = Ynj*a;
g = reshape(g, m, n);

% Visualise fit.
subplot(1, 3, 2);
surf(x, y, z, g, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
colorbar;
daspect([1, 1, 1]);
view(3);

% Visualise difference.
subplot(1, 3, 3);
surf(x, y, z, f-g, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
colorbar;
daspect([1, 1, 1]);
view(3);

fprintf('Least-squares error is %f.\n', norm(f-g, 2)^2);

end