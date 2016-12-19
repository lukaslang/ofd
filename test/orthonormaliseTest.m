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
function tests = orthonormaliseTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create two test vectors.
v = [3, 1, 0];
w = [2, 2, 0];

% Orthonormalise.
[A, e, f] = orthonormalise(v, w);

% Check first basis vector.
verifyEqual(testCase, e, [3/sqrt(10), 1/sqrt(10), 0], 'AbsTol', 1e-15);
% Check second basis vector.
verifyEqual(testCase, f, [-1/sqrt(10), 3/sqrt(10), 0], 'AbsTol', 1e-15);

% Check lengths.
verifyEqual(testCase, sqrt(sum(e .^ 2, 2)), 1, 'AbsTol', 1e-15);
verifyEqual(testCase, sqrt(sum(f .^ 2, 2)), 1, 'AbsTol', 1e-15);

% Check orthogonality.
verifyEqual(testCase, e * f', 0, 'AbsTol', 1e-15);

% Check coefficients.
verifyEqual(testCase, A(1, 1, :), 1/sqrt(10), 'AbsTol', 1e-15);
verifyEqual(testCase, A(1, 2, :), 0, 'AbsTol', 1e-15);
verifyEqual(testCase, A(2, 1, :), -8/(10 * sqrt(4/25 + 36/25)), 'AbsTol', 1e-15);
verifyEqual(testCase, A(2, 2, :), 1/sqrt(4/25 + 36/25), 'AbsTol', 1e-15);

% Check that A transforms {v, w} into {e, f}.
verifyEqual(testCase, A * [v; w], [e; f], 'AbsTol', 1e-15);

end

function perturbedSphereTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0:2;
Y = spharm(0, [0, 0, 1]);
c = [1 / Y, 0, 0, 0, 0, 0, 0, 0.5, 0]';

% Compute synthesis.
S = surfsynth(Ns, V, c);

% Compute tangential basis at incenters.
[Z, IC] = surftangentialbasis(Ns, c, F, V);
Z1 = squeeze(Z(:, 1, :));
Z2 = squeeze(Z(:, 2, :));

% Compute orthonormal basis.
[~, e, f] = orthonormalise(Z1, Z2);

% Check if orthogonal.
ip = sum(e .* f, 2);
verifyEqual(testCase, zeros(m, 1), ip, 'AbsTol', 1e-12);

% Check length one.
verifyEqual(testCase, sqrt(sum(e.^2, 2)), ones(m, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, sqrt(sum(f.^2, 2)), ones(m, 1), 'AbsTol', 1e-12);

% Plot tangential basis and dot product.
figure;
axis([-1, 1, -1, 1, -1, 1]);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), sum(Z1 .* Z2, 2));
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z1(:, 1), Z1(:, 2), Z1(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), Z2(:, 1), Z2(:, 2), Z2(:, 3), 1, 'b');

figure;
axis([-1, 1, -1, 1, -1, 1]);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), ip);
daspect([1, 1, 1]);
colorbar;
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), e(:, 1), e(:, 2), e(:, 3), 1, 'r');
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), f(:, 1), f(:, 2), f(:, 3), 1, 'b');
end