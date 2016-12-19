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
function tests = triangsectionTest
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

[Fs, Vs] = triangsection(F, V, [-1, 1, -1, 1, 0, 1]);

verifyEqual(testCase, Fs, F);
verifyEqual(testCase, Vs, V);

end

function visualiseTest(testCase)

% Generate icosahedron.
tic;
[F, V] = halfsphTriang(4);
toc;
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
assert(all(V(:, 3) >= 0));

[Fs, Vs] = triangsection(F, V, [-1, 0, 0, 1, 0, 1]);

figure;
trisurf(Fs, Vs(:, 1), Vs(:, 2), Vs(:, 3));
daspect([1, 1, 1]);

end