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
function test_suite = vspharmTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

% Create spherical harmonics.
N = 5;
[Y1, Y2] = vspharm(N, F, V);
assertFalse(isempty(Y1));
assertEqual(size(Y1), [n, 2*N + 1, 3]);
assertFalse(isempty(Y2));
assertEqual(size(Y2), [n, 2*N + 1, 3]);

end

function orthogonalityTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

% Create spherical harmonics.
N = 5;
[Y1, Y2] = vspharm(N, F, V);

% Compute R3 inner product.
ip = dot(Y1, Y2, 3);
assertAlmostEqual(ip, zeros(n, 2*N + 1));

end

function visualiseTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

T = TriRep(F, V);
P = T.incenters;

% Create spherical harmonics.
N = 3;
[Y1, Y2] = vspharm(N, F, V);
assertFalse(isempty(Y1));
assertEqual(size(Y1), [n, 2*N + 1, 3]);
assertFalse(isempty(Y2));
assertEqual(size(Y2), [n, 2*N + 1, 3]);

% Create spherical harmonics for visualisation.
Ynj = spharm(N, V);

figure;
for k=1:2*N+1
    f = Ynj(:, k);
    subplot(1, 2*N+1, k);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(P(:, 1), P(:, 2), P(:, 3), Y1(:, k, 1), Y1(:, k, 2), Y1(:, k, 3), 1, 'r');
    quiver3(P(:, 1), P(:, 2), P(:, 3), Y2(:, k, 1), Y2(:, k, 2), Y2(:, k, 3), 1, 'b');
end
end