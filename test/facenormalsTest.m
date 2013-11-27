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
function test_suite = facenormalsTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
V(1, :) = [0, 0, 0];
V(2, :) = [0, 2, 0];
V(3, :) = [3, 0, 0];
% Orient triangle clockwise.
F = 1:3;

% Compute face unit normals.
Fn = facenormals(F, V);
assertFalse(isempty(Fn));
assertEqual(size(Fn), [1, 3]);
assertAlmostEqual(Fn, [0, 0, 1]);

% Check length.
assertAlmostEqual(sqrt(sum(Fn.^2, 2)), 1);

end

function sphTriangTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Compute face unit normals.
Fn = facenormals(F, V);
assertFalse(isempty(Fn));
assertEqual(size(Fn), [n, 3]);

% Check length.
assertAlmostEqual(sqrt(sum(Fn.^2, 2)), ones(n, 1));

end

function visualisationTest

% Generate icosahedron.
[F, V] = sphTriang(2);
n = size(F, 1);
assertFalse(isempty(F));
assertFalse(isempty(V));

% Compute face unit normals.
Fn = facenormals(F, V);
assertFalse(isempty(Fn));

% Check length.
assertAlmostEqual(sqrt(sum(Fn.^2, 2)), ones(n, 1));

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));

% Get incenters of triangles.
TR = TriRep(F, V);
P = TR.incenters;

% Plot vectors.
scale = 0;
quiver3(P(:, 1), P(:, 2), P(:, 3), Fn(:, 1), Fn(:, 2), Fn(:, 3), scale);
daspect([1, 1, 1]);
view(3);

end