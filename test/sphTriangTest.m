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
function test_suite = sphTriangTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = sphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [20, 3]);
assertEqual(size(V), [12, 3]);

% Generate icosahedron.
[F, V] = sphTriang(0);
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [20, 3]);
assertEqual(size(V), [12, 3]);

end

function unitSphereTest

% Generate icosahedron.
[~, V] = sphTriang;

% Check if vertices lie on unit sphere.
assertAlmostEqual(sqrt(sum(V.^2, 2)), ones(12, 1));

end

function resultWithRefinementTest

% Generate icosahedron.
[F, V] = sphTriang(1);
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [80, 3]);
assertEqual(size(V), [42, 3]);

% Check if vertices lie on unit sphere.
assertAlmostEqual(sqrt(sum(V.^2, 2)), ones(42, 1));

% Generate icosahedron.
[F, V] = sphTriang(3);
assertFalse(isempty(F));
assertFalse(isempty(V));

% Check if vertices lie on unit sphere.
assertAlmostEqual(sqrt(sum(V.^2, 2)), ones(size(V, 1), 1));

end

function createTriRepTest

% Generate icosahedron.
[F, V] = sphTriang;
T = TriRep(F, V);
assertFalse(isempty(T));
assertEqual(T.size, [20, 3]);

end

function visualiseTest

% Generate icosahedron.
[F, V] = sphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [20, 3]);
assertEqual(size(V), [12, 3]);

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
daspect([1, 1, 1]);

end

function visualiseRefinementTest

% Generate icosahedron.
tic;
[F, V] = sphTriang(4);
toc;
assertFalse(isempty(F));
assertFalse(isempty(V));

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
daspect([1, 1, 1]);

end