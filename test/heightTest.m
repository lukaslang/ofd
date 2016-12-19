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
function tests = heightTest
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
[F, V] = sphTriang;
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [20, 3]);
verifyEqual(testCase, size(V), [12, 3]);

% Compute height vectors.
H = height(F, V);
verifyFalse(testCase, isempty(H));
verifyEqual(testCase, size(H), [20, 3, 3]);

end

function visualisationTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang(3);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));

% Compute height vectors.
H = height(F, V);
verifyFalse(testCase, isempty(H));
verifyEqual(testCase, size(H), [size(F, 1), 3, 3]);

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
% Get coordinates of vertices.
V1 = V(F(:, 1), :);
V2 = V(F(:, 2), :);
V3 = V(F(:, 3), :);
% Plot vectors.
scale = 1;
quiver3(V1(:, 1), V1(:, 2), V1(:, 3), H(:, 1, 1), H(:, 1, 2), H(:, 1, 3), scale);
quiver3(V2(:, 1), V2(:, 2), V2(:, 3), H(:, 2, 1), H(:, 2, 2), H(:, 2, 3), scale);
quiver3(V3(:, 1), V3(:, 2), V3(:, 3), H(:, 3, 1), H(:, 3, 2), H(:, 3, 3), scale);
daspect([1, 1, 1]);
view(3);

end