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
function tests = surfpushforwardTest
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
[F, V] = sphTriang(3);
m = size(F, 1);
n = size(V, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Compute synthesis.
[Vs, rho] = surfsynth(Ns, V, c);

% Check if surface is unit sphere.
len = sqrt(sum(Vs .^ 2, 2));
verifyEqual(testCase, len, ones(n, 1), 'AbsTol', 1e-12);

% Check if function rho ise identically one.
verifyEqual(testCase, rho, ones(n, 1), 'AbsTol', 1e-12);

% Compute gradient of rho on triangulation.
g = grad(F, V, rho);
% Check if gradient is zero.
verifyEqual(testCase, g, zeros(m, 3), 'AbsTol', 1e-12);

% Compute incenters.
T = TriRep(F, V);
IC = T.incenters;

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, Y2] = vspharm(deg, F, V);

% Compute the pushforward.
[Z1, ICs] = surfpushforward(Ns, c, F, V, Y1);
verifyEqual(testCase, Y1, Z1, 'AbsTol', 1e-12);
verifyEqual(testCase, IC, ICs, 'AbsTol', 5e-3);

% Compute the pushforward.
[Z2, ICs] = surfpushforward(Ns, c, F, V, Y2);
verifyEqual(testCase, Y2, Z2, 'AbsTol', 1e-12);
verifyEqual(testCase, IC, ICs, 'AbsTol', 5e-3);

end

function identityMapTangentialTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
m = size(F, 1);

% Compute face normals.
T = TriRep(F, V);
FN = -T.faceNormals;
IC = T.incenters;

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, Y2] = vspharm(deg, F, V);

% Compute the pushforward.
[Z1, ICs] = surfpushforward(Ns, c, F, V, Y1);
verifyEqual(testCase, Y1, Z1, 'AbsTol', 1e-12);
verifyEqual(testCase, IC, ICs, 'AbsTol', 3e-4);

% Compute the pushforward.
[Z2, ICs] = surfpushforward(Ns, c, F, V, Y2);
verifyEqual(testCase, Y2, Z2, 'AbsTol', 1e-12);
verifyEqual(testCase, IC, ICs, 'AbsTol', 3e-4);

% Check if orthogonal to surface normal.
for k=1:size(Z1, 2)
    verifyEqual(testCase, dot(squeeze(Z1(:, k, :)), FN, 2), zeros(m, 1), 'AbsTol', 1e-12);
    verifyEqual(testCase, dot(squeeze(Z2(:, k, :)), FN, 2), zeros(m, 1), 'AbsTol', 1e-12);
end
end

function scaledSphereTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(V, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 2 / Y;

% Compute synthesis.
[Vs, rho] = surfsynth(Ns, V, c);

% Check if surface is sphere with radius 2.
len = sqrt(sum(Vs .^ 2, 2));
verifyEqual(testCase, len, 2*ones(n, 1), 'AbsTol', 1e-12);

% Check if function rho ise identically two.
verifyEqual(testCase, rho, 2*ones(n, 1), 'AbsTol', 1e-12);

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, Y2] = vspharm(deg, F, V);

% Compute the pushforwards.
[Z1, ~] = surfpushforward(Ns, c, F, V, Y1);
[Z2, ICs] = surfpushforward(Ns, c, F, V, Y2);
verifyEqual(testCase, Z1, 2 * Y1, 'AbsTol', 1e-12);
verifyEqual(testCase, Z2, 2 * Y2, 'AbsTol', 1e-12);

% Create spherical harmonics for visualisation.
Ynj = spharm(deg, V);

figure;
for k=1:2*deg+1
    f = Ynj(:, k);
    subplot(1, 2*deg+1, k);
    axis([-2, 2, -2, 2, -2, 2]);
    hold on;
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(ICs(:, 1), ICs(:, 2), ICs(:, 3), Z1(:, k, 1), Z1(:, k, 2), Z1(:, k, 3), 1, 'r');
    quiver3(ICs(:, 1), ICs(:, 2), ICs(:, 3), Z2(:, k, 1), Z2(:, k, 2), Z2(:, k, 3), 1, 'b');
end

end

function scaledSphereTangentialTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);
n = size(V, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 2 / Y;

% Compute synthesis.
[Vs, rho] = surfsynth(Ns, V, c);

% Check if surface is sphere with radius 2.
len = sqrt(sum(Vs .^ 2, 2));
verifyEqual(testCase, len, 2*ones(n, 1), 'AbsTol', 1e-12);

% Check if function rho ise identically two.
verifyEqual(testCase, rho, 2*ones(n, 1), 'AbsTol', 1e-12);

% Compute face normals.
T = TriRep(F, Vs);
FN = -T.faceNormals;

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, Y2] = vspharm(deg, F, V);

% Compute the pushforwards.
[Z1, ~] = surfpushforward(Ns, c, F, V, Y1);
[Z2, ~] = surfpushforward(Ns, c, F, V, Y2);
verifyEqual(testCase, Z1, 2 * Y1, 'AbsTol', 1e-12);
verifyEqual(testCase, Z2, 2 * Y2, 'AbsTol', 1e-12);

% Check if orthogonal to surface normal.
for k=1:size(Z1, 2)
    verifyEqual(testCase, dot(squeeze(Z1(:, k, :)), FN, 2), zeros(m, 1), 'AbsTol', 1e-12);
    verifyEqual(testCase, dot(squeeze(Z2(:, k, :)), FN, 2), zeros(m, 1), 'AbsTol', 1e-12);
end

end

function perturbedSphereTangentialTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
m = size(F, 1);

% Set parameters for perturbed sphere.
Ns = 0:2;
Y = spharm(0, [0, 0, 1]);
c = [1 / Y, 0, 0, 0, 0, 0, 0, 0.5, 0]';

% Compute synthesis.
[Vs, ~] = surfsynth(Ns, V, c);

% Compute face normals.
T = TriRep(F, Vs);
FN = -T.faceNormals;

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, Y2] = vspharm(deg, F, V);

% Compute the pushforwards.
[Z1, ~] = surfpushforward(Ns, c, F, V, Y1);
[Z2, ~] = surfpushforward(Ns, c, F, V, Y2);

% Check if orthogonal to surface normal.
for k=1:size(Z1, 2)
    verifyEqual(testCase, dot(squeeze(Z1(:, k, :)), FN, 2), zeros(m, 1), 'AbsTol', 2e-3);
    verifyEqual(testCase, dot(squeeze(Z2(:, k, :)), FN, 2), zeros(m, 1), 'AbsTol', 2e-3);
end

end

function perturbedSphereVisualisationTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);

% Set parameters for perturbed sphere.
Ns = 0:2;
Y = spharm(0, [0, 0, 1]);
c = [1 / Y, 0, 0, 0, 0, 0, 0, 0.5, 0]';

% Compute synthesis.
[Vs, ~] = surfsynth(Ns, V, c);

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, Y2] = vspharm(deg, F, V);

% Compute the pushforwards.
[Z1, ~] = surfpushforward(Ns, c, F, V, Y1);
[Z2, ICs] = surfpushforward(Ns, c, F, V, Y2);

% Create spherical harmonics for visualisation.
Ynj = spharm(deg, V);

figure;
for k=1:2*deg+1
    f = Ynj(:, k);
    subplot(1, 2*deg+1, k);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(ICs(:, 1), ICs(:, 2), ICs(:, 3), Z1(:, k, 1), Z1(:, k, 2), Z1(:, k, 3), 1, 'r');
    quiver3(ICs(:, 1), ICs(:, 2), ICs(:, 3), Z2(:, k, 1), Z2(:, k, 2), Z2(:, k, 3), 1, 'b');
end

end