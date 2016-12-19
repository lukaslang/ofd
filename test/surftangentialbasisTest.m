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
function tests = surftangentialbasisTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function identityMapTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
n = size(V, 1);
m = size(F, 1);

% Set parameters for unit sphere.
N = 0;
Y = spharm(N, [0, 0, 1]);
c = 1 / Y;

% Compute synthesis.
[Vs, rho] = surfsynth(N, V, c);

% Check if surface is unit sphere.
len = sqrt(sum(Vs .^ 2, 2));
verifyEqual(testCase, len, ones(n, 1), 'AbsTol', 1e-12);

% Check if function rho is identically one.
verifyEqual(testCase, rho, ones(n, 1), 'AbsTol', 1e-12);

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(N, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Check if orthogonal.
ip = sum(Z1 .* Z2, 2);
verifyEqual(testCase, zeros(m, 1), ip, 'AbsTol', 1e-12);

% Check length one.
verifyEqual(testCase, sqrt(sum(Z1.^2, 2)), ones(m, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, sqrt(sum(Z2.^2, 2)), ones(m, 1), 'AbsTol', 1e-12);

% Plot tangential basis and dot product.
figure;
axis([-1, 1, -1, 1, -1, 1]);
hold on;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), ip);
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');
end

function identityMapIncentersTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);

% Set parameters for unit sphere.
N = 0;
Y = spharm(N, [0, 0, 1]);
c = 1 / Y;

% Compute tangential basis at incenters.
[~, IC] = surftangentialbasis(N, c, F, V);

% Check if orthogonal to surface normal.
T = TriRep(F, V);
len = sqrt(sum(T.incenters .^2, 2));
verifyEqual(testCase, bsxfun(@rdivide, T.incenters, len), IC, 'AbsTol', 1e-12);

end

function identityMapTangentialTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
m = size(F, 1);

% Compute face normals.
T = TriRep(F, V);
FN = -T.faceNormals;

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(Ns, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Plot tangential basis and dot product with face normals.
figure;
subplot(1, 2, 1);
axis([-1, 1, -1, 1, -1, 1]);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dot(Z1, FN, 2));
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), FN(:, 1), FN(:, 2), FN(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'b');

subplot(1, 2, 2);
axis([-1, 1, -1, 1, -1, 1]);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dot(Z2, FN, 2));
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), FN(:, 1), FN(:, 2), FN(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');

% Check if orthogonal to surface normal.
warning('Note that surftangentialbasis does not return vectors tangent to face normals!');
verifyEqual(testCase, dot(Z1, FN, 2), zeros(m, 1), 'AbsTol', 4e-3);
verifyEqual(testCase, dot(Z2, FN, 2), zeros(m, 1), 'AbsTol', 4e-3);

end

function scaledSphereTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
n = size(V, 1);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 2 / Y;

% Compute synthesis.
[Vs, rho] = surfsynth(Ns, V, c);

% Check if surface is sphere of radius two.
len = sqrt(sum(Vs .^ 2, 2));
verifyEqual(testCase, len, 2 * ones(n, 1), 'AbsTol', 1e-12);

% Check if function rho ise identically two.
verifyEqual(testCase, rho, 2 * ones(n, 1), 'AbsTol', 1e-12);

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(Ns, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Check if orthogonal.
ip = sum(Z1 .* Z2, 2);
verifyEqual(testCase, ip, zeros(m, 1), 'AbsTol', 1e-12);

% Plot tangential basis and dot product.
figure;
axis([-2, 2, -2, 2, -2, 2]);
hold on;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), ip);
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');
end

function perturbedSphereVisualisationTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0:2;
Y = spharm(0, [0, 0, 1]);
c = [1 / Y, 0, 0, 0, 0, 0, 0, 0.5, 0]';

% Compute synthesis.
Vs = surfsynth(Ns, V, c);

% Compute face normals.
T = TriRep(F, Vs);
FN = -T.faceNormals;

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(Ns, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Compute dot product.
ip = sum(Z1 .* Z2, 2);

% Plot tangential basis and dot product.
figure;
axis([-2, 2, -2, 2, -2, 2]);
hold on;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), ip);
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');

% Check if orthogonal to surface normal.
warning('Note that surftangentialbasis does not return vectors tangent to face normals!');
verifyEqual(testCase, dot(Z1, FN, 2), zeros(m, 1), 'AbsTol', 7e-3);
verifyEqual(testCase, dot(Z2, FN, 2), zeros(m, 1), 'AbsTol', 7e-3);

end