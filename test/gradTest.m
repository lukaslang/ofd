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
function test_suite = gradTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = sphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [20, 3]);
assertEqual(size(V), [12, 3]);

% Create sample function.
f = ones(12, 1);

% Compute gradient.
g = grad(F, V, f);
assertFalse(isempty(g));
assertEqual(size(g), [20, 3]);
assertEqual(g, zeros(20, 3));

end

function visualiseTest

% Generate icosahedron.
[F, V] = sphTriang(3);
assertFalse(isempty(F));
assertFalse(isempty(V));

T = TriRep(F, V);
P = T.incenters;

% Create spherical harmonics as test functions.
N = 3;
f = spharm(N, V);

figure;
for k=1:2*N+1
    % Compute gradient of k-th harmonic.
    g = grad(F, V, f(:, k));
    assertFalse(isempty(g));

    % Plot spherical harmonics.
    subplot(1, 2*N+1, k);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f(:, k));
    shading interp;
    daspect([1, 1, 1]);
    view(3);

    % Plot gradient field.
    quiver3(P(:, 1), P(:, 2), P(:, 3), g(:, 1), g(:, 2), g(:, 3));
end
end