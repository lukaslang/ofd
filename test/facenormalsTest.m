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
function tests = facenormalsTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create triangulation of unit sphere.
V(1, :) = [0, 0, 0];
V(2, :) = [0, 2, 0];
V(3, :) = [3, 0, 0];
% Orient triangle clockwise.
F = 1:3;

% Compute face unit normals.
Fn = facenormals(F, V);
verifyFalse(testCase, isempty(Fn));
verifyEqual(testCase, size(Fn), [1, 3]);
verifyEqual(testCase, Fn, [0, 0, 1], 'AbsTol', 1e-15);

% Check length.
verifyEqual(testCase, sqrt(sum(Fn.^2, 2)), 1, 'AbsTol', 1e-15);

end

function sphTriangTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Compute face unit normals.
Fn = facenormals(F, V);
verifyFalse(testCase, isempty(Fn));
verifyEqual(testCase, size(Fn), [n, 3]);

% Check length.
verifyEqual(testCase, sqrt(sum(Fn.^2, 2)), ones(n, 1), 'AbsTol', 1e-15);

end

function visualisationTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang(2);
n = size(F, 1);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));

% Compute face unit normals.
Fn = facenormals(F, V);
verifyFalse(testCase, isempty(Fn));

% Check length.
verifyEqual(testCase, sqrt(sum(Fn.^2, 2)), ones(n, 1), 'AbsTol', 1e-15);

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