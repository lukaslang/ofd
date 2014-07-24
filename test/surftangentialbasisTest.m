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
function test_suite = surftangentialbasisTest
    initTestSuite;
end

function identityMapTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
n = size(V, 1);
m = size(F, 1);

% Set parameters for unit sphere.
N = 0;
Y = spharm(N, [0, 0, 1]);
c = 1 / Y;

% Compute synthesis.
[S, rho] = surfsynth(N, V, c);

% Check if surface is unit sphere.
len = sqrt(sum(S .^ 2, 2));
assertAlmostEqual(len, ones(n, 1), 1e-12);

% Check if function rho ise identically one.
assertAlmostEqual(rho, ones(n, 1), 1e-12);

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(N, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Check if orthogonal.
ip = sum(Z1 .* Z2, 2);
assertAlmostEqual(zeros(m, 1), ip, 1e-12);

% Check length one.
assertAlmostEqual(sqrt(sum(Z1.^2, 2)), ones(m, 1), 1e-12);
assertAlmostEqual(sqrt(sum(Z2.^2, 2)), ones(m, 1), 1e-12);

% Plot tangential basis and dot product.
figure;
axis([-1, 1, -1, 1, -1, 1]);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), ip);
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');
end

function scaledSphereTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
n = size(V, 1);
m = size(F, 1);

% Set parameters for unit sphere.
N = 0;
Y = spharm(N, [0, 0, 1]);
c = 2 / Y;

% Compute synthesis.
[S, rho] = surfsynth(N, V, c);

% Check if surface is sphere of radius two.
len = sqrt(sum(S .^ 2, 2));
assertAlmostEqual(len, 2 * ones(n, 1), 1e-12);

% Check if function rho ise identically two.
assertAlmostEqual(rho, 2 * ones(n, 1), 1e-12);

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(N, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Check if orthogonal.
ip = sum(Z1 .* Z2, 2);
assertAlmostEqual(zeros(m, 1), ip, 1e-12);

% Plot tangential basis and dot product.
figure;
axis([-2, 2, -2, 2, -2, 2]);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), ip);
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');
end

function perturbedSphereVisualisationTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);

% Set parameters for unit sphere.
N = 0:2;
Y = spharm(0, [0, 0, 1]);
c = [1 / Y, 0, 0, 0, 0, 0, 0, 0.5, 0]';

% Compute synthesis.
S = surfsynth(N, V, c);

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(N, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Compute dot product.
ip = sum(Z1 .* Z2, 2);

% Plot tangential basis and dot product.
figure;
axis([-2, 2, -2, 2, -2, 2]);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), ip);
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');
end