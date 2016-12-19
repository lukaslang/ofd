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
function tests = refineMeshTest
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

numF = size(F, 1);
numV = size(V, 1);

% Refine.
[F, V] = refineMesh(F, V);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [80, 3]);
verifyEqual(testCase, size(V), [42, 3]);

end

function visualiseTest(testCase)

[F, V] = sphTriang;
[F, V] = refineMesh(F, V);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [80, 3]);
verifyEqual(testCase, size(V), [42, 3]);

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3));
daspect([1, 1, 1]);

end