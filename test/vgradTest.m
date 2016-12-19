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
function tests = vgradTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function vgradZeroVectorFieldTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
m = size(F, 1);
n = size(V, 1);

% Create a zero piecewise constant vector field on the unit sphere.
Y = zeros(n, 3);

% Compute vectorial gradient of Y.
G = vgrad(F, V, Y, ones(m, 3));
verifyEqual(testCase, G, zeros(m, 3), 'AbsTol', 1e-12);

end

function vgradHeightTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
m = size(F, 1);
n = size(V, 1);

% Create a zero piecewise constant vector field on the unit sphere.
Y = zeros(n, 3);

H = height(F, V);

% Compute vectorial gradient of Y.
G = vgrad(F, V, Y, ones(m, 3), H);
verifyEqual(testCase, G, zeros(m, 3), 'AbsTol', 1e-12);

end

function vgradlenHTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
m = size(F, 1);
n = size(V, 1);

% Create a zero piecewise constant vector field on the unit sphere.
Y = zeros(n, 3);

H = height(F, V);
lenH = sum(H(:, 2:3, :).^2, 3);

% Compute vectorial gradient of Y.
G = vgrad(F, V, Y, ones(m, 3), H, lenH);
verifyEqual(testCase, G, zeros(m, 3), 'AbsTol', 1e-12);

end

function vgradFaceNormalsTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
m = size(F, 1);
n = size(V, 1);

% Create a zero piecewise constant vector field on the unit sphere.
Y = zeros(n, 3);

H = height(F, V);
lenH = sum(H(:, 2:3, :).^2, 3);
T = TriRep(F, V);
FN = -T.faceNormals;

% Compute vectorial gradient of Y.
G = vgrad(F, V, Y, ones(m, 3), H, lenH, FN);
verifyEqual(testCase, G, zeros(m, 3), 'AbsTol', 1e-12);

end

% TODO: Write testcase for tangent vgrad!
function vgradTangentTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
m = size(F, 1);

% Compute triangulation of incenters.
T = TriRep(F, V);
DV = T.incenters;
% Compute triangles spanned by incenters of neighboring faces.
DF = T.neighbors;

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;
% Compute tangential basis at incenters.
[Z, ~] = surftangentialbasis(Ns, c, F, V);
Z = squeeze(Z(:, 1, :));

% Project to tangent space.
FN = -T.faceNormals;
Z = Z - bsxfun(@times, FN, dot(Z, FN, 2));
verifyEqual(testCase, dot(Z, FN, 2), zeros(m, 1), 'AbsTol', 1e-15);

% Create piecewise constant vector field on the unit sphere.
deg = 2;
k = 1;
[Y, ~] = vspharm(deg, F, V);
Y = squeeze(Y(:, k, :));

% Compute covariant derivative.
G = vgrad(DF, DV, Y, Z);

% Check if tangential to given triangulation (DF, DV).
DT = TriRep(DF, DV);
FN = -DT.faceNormals;
verifyEqual(testCase, dot(G, FN, 2), zeros(m, 1), 'AbsTol', 1e-15);

% Check if tangential to original triangulation (F, V).
FN = -T.faceNormals;
verifyEqual(testCase, dot(G, FN, 2), zeros(m, 1), 'AbsTol', 7e-3);

end

function vgradSphereVisualisationTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, ~] = vspharm(deg, F, V);

% Compute triangulation of incenters.
T = TriRep(F, V);
DV = T.incenters;
% Compute triangles spanned by incenters of neighboring faces.
DF = T.neighbors;

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;
% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(Ns, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));
% Compute orthonormal basis.
[~, Z1, Z2] = orthonormalise(Z1, Z2);

% Create spherical harmonics for visualisation.
Ynj = spharm(deg, V);

for k=1:2*deg+1
    % Visualise type 2 vector spherical harmonics.
    Y = squeeze(Y1(:, k, :));
    % Compute gradient of Y.
    G = vgrad(DF, DV, Y, Z1);
    f = Ynj(:, k);
    figure;
    subplot(1, 2, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Y(:, 1), Y(:, 2), Y(:, 3), 1, 'r');
    % Plot second vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'y');
    % Plot gradient.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), G(:, 1), G(:, 2), G(:, 3), 1, 'b');
    subplot(1, 2, 2);
    % Compute squared norm of gradient.
    f = sum(G .^2, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    daspect([1, 1, 1]);
    view(3);
end

for k=1:2*deg+1
    % Visualise type 2 vector spherical harmonics.
    Y = squeeze(Y1(:, k, :));
    % Compute gradient of Y.
    G = vgrad(DF, DV, Y, Z2);
    f = Ynj(:, k);
    figure;
    subplot(1, 2, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Y(:, 1), Y(:, 2), Y(:, 3), 1, 'r');
    % Plot second vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'y');
    % Plot gradient.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), G(:, 1), G(:, 2), G(:, 3), 1, 'b');
    subplot(1, 2, 2);
    % Compute squared norm of gradient.
    f = sum(G .^2, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    daspect([1, 1, 1]);
    view(3);
end

end

function perturbedSphereVisualisationTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
m = size(F, 1);

% Create piecewise constant vector field on the unit sphere.
deg = 2;
[Y1, ~] = vspharm(deg, F, V);

% Set parameters for perturbed sphere.
Ns = 0:2;
Y = spharm(0, [0, 0, 1]);
c = [1 / Y, 0, 0, 0, 0, 0, 0, 0.5, 0]';

% Compute synthesis.
Vs = surfsynth(Ns, V, c);

% Compute triangulation of incenters.
T = TriRep(F, Vs);
FN = -T.faceNormals;
DV = T.incenters;
% Compute triangles spanned by incenters of neighboring faces.
DF = T.neighbors;

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(Ns, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));
verifyEqual(testCase, dot(Z1, FN, 2), zeros(m, 1), 'AbsTol', 7e-3);
verifyEqual(testCase, dot(Z2, FN, 2), zeros(m, 1), 'AbsTol', 7e-3);

% Compute orthonormal basis.
[~, E1, E2] = orthonormalise(Z1, Z2);
verifyEqual(testCase, dot(E1, FN, 2), zeros(m, 1), 'AbsTol', 7e-3);
verifyEqual(testCase, dot(E2, FN, 2), zeros(m, 1), 'AbsTol', 7e-3);

% Check if orthogonal.
ip = sum(E1 .* E2, 2);
verifyEqual(testCase, zeros(m, 1), ip, 'AbsTol', 1e-12);

% Create spherical harmonics for visualisation.
Ynj = spharm(deg, V);

for k=1:2*deg+1
    % Visualise type 2 vector spherical harmonics.
    Y = squeeze(Y1(:, k, :));
    % Compute gradient of Y.
    G = vgrad(DF, DV, Y, E1);    
    f = Ynj(:, k);
    figure;
    subplot(1, 2, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Y(:, 1), Y(:, 2), Y(:, 3), 1, 'r');
    % Plot second vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), E1(:, 1), E1(:, 2), E1(:, 3), 1, 'y');
    % Plot gradient.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), G(:, 1), G(:, 2), G(:, 3), 1, 'b');
    subplot(1, 2, 2);
    % Compute squared norm of gradient.
    f = sum(G .^2, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f);
    daspect([1, 1, 1]);
    view(3);
end

for k=1:2*deg+1
    % Visualise type 2 vector spherical harmonics.
    Y = squeeze(Y1(:, k, :));
    % Compute gradient of Y.
    G = vgrad(DF, DV, Y, E2);
    f = Ynj(:, k);
    figure;
    subplot(1, 2, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Y(:, 1), Y(:, 2), Y(:, 3), 1, 'r');
    % Plot second vector field.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), E2(:, 1), E2(:, 2), E2(:, 3), 1, 'y');
    % Plot gradient.
    quiver3(IC(:, 1), IC(:, 2), IC(:, 3), G(:, 1), G(:, 2), G(:, 3), 1, 'b');
    subplot(1, 2, 2);
    % Compute squared norm of gradient.
    f = sum(G .^2, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f);
    daspect([1, 1, 1]);
    view(3);
end

end