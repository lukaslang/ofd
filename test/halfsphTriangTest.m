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
function tests = halfsphTriangTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Generate icosahedron.
[F, V] = halfsphTriang;
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [8, 3]);
verifyEqual(testCase, size(V), [8, 3]);
assert(all(V(:, 3) >= 0));

% Generate icosahedron.
[F, V] = halfsphTriang(0);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [8, 3]);
verifyEqual(testCase, size(V), [8, 3]);
assert(all(V(:, 3) >= 0));

end

function resultWithRefinementTest(testCase)

% Generate icosahedron.
[F, V] = halfsphTriang(1);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [36, 3]);
verifyEqual(testCase, size(V), [25, 3]);
assert(all(V(:, 3) >= 0));

% Check if vertices lie on unit sphere.
verifyEqual(testCase, sqrt(sum(V.^2, 2)), ones(25, 1), 'AbsTol', 1e-15);

% Generate icosahedron.
[F, V] = halfsphTriang(3);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
assert(all(V(:, 3) >= 0));

% Check if vertices lie on unit sphere.
verifyEqual(testCase, sqrt(sum(V.^2, 2)), ones(size(V, 1), 1), 'AbsTol', 1e-15);

end

function createTriRepTest(testCase)

% Generate icosahedron.
[F, V] = halfsphTriang;
T = TriRep(F, V);
verifyFalse(testCase, isempty(T));
verifyEqual(testCase, T.size, [8, 3]);

end

function visualiseTest(testCase)

% Generate icosahedron.
[F, V] = halfsphTriang;
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [8, 3]);
verifyEqual(testCase, size(V), [8, 3]);
assert(all(V(:, 3) >= 0));

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
daspect([1, 1, 1]);

end

function visualiseRefinementTest(testCase)

% Generate icosahedron.
tic;
[F, V] = halfsphTriang(4);
toc;
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
assert(all(V(:, 3) >= 0));

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
daspect([1, 1, 1]);

end