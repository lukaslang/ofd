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
function test_suite = halfsphTriangTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = halfsphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [8, 3]);
assertEqual(size(V), [8, 3]);
assert(all(V(:, 3) >= 0));

% Generate icosahedron.
[F, V] = halfsphTriang(0);
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [8, 3]);
assertEqual(size(V), [8, 3]);
assert(all(V(:, 3) >= 0));

end

function resultWithRefinementTest

% Generate icosahedron.
[F, V] = halfsphTriang(1);
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [36, 3]);
assertEqual(size(V), [25, 3]);
assert(all(V(:, 3) >= 0));

% Check if vertices lie on unit sphere.
assertAlmostEqual(sqrt(sum(V.^2, 2)), ones(25, 1));

% Generate icosahedron.
[F, V] = halfsphTriang(3);
assertFalse(isempty(F));
assertFalse(isempty(V));
assert(all(V(:, 3) >= 0));

% Check if vertices lie on unit sphere.
assertAlmostEqual(sqrt(sum(V.^2, 2)), ones(size(V, 1), 1));

end

function createTriRepTest

% Generate icosahedron.
[F, V] = halfsphTriang;
T = TriRep(F, V);
assertFalse(isempty(T));
assertEqual(T.size, [8, 3]);

end

function visualiseTest

% Generate icosahedron.
[F, V] = halfsphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [8, 3]);
assertEqual(size(V), [8, 3]);
assert(all(V(:, 3) >= 0));

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
daspect([1, 1, 1]);

end

function visualiseRefinementTest

% Generate icosahedron.
tic;
[F, V] = halfsphTriang(7);
toc;
assertFalse(isempty(F));
assertFalse(isempty(V));
assert(all(V(:, 3) >= 0));

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
daspect([1, 1, 1]);

end