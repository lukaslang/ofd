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
function tests = linearsystemTest
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
[Faces, Verts] = sphTriang(3);
m = size(Verts, 1);

% Create random data.
f1 = randi(255, m, 1);
f2 = randi(255, m, 1);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
N = 10;
[dim, U, d, b] = linearsystem(Faces, Verts, 1:N, f1, f2, h, tol);

% Check results.
verifyEqual(testCase, dim, 2*(N^2 + 2*N));
verifyFalse(testCase, isempty(U));
verifyEqual(testCase, size(U), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(d));
verifyTrue(testCase, isvector(d));
verifyEqual(testCase, size(d), [2*(N^2 + 2*N), 1]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [2*(N^2 + 2*N), 1]);

% Check if matrix is symmetric.
verifyEqual(testCase, U, U', 'AbsTol', 1e-15);

end

function intervalTest(testCase)

% Create triangulation of unit sphere.
[Faces, Verts] = sphTriang(3);
m = size(Verts, 1);

% Create random data.
f1 = randi(255, m, 1);
f2 = randi(255, m, 1);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
N = 3:10;
[dim, U, d, b] = linearsystem(Faces, Verts, N, f1, f2, h, tol);

% Check results.
expDim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);
verifyEqual(testCase, dim, expDim);
verifyFalse(testCase, isempty(U));
verifyEqual(testCase, size(U), [expDim, expDim]);
verifyFalse(testCase, isempty(d));
verifyTrue(testCase, isvector(d));
verifyEqual(testCase, size(d), [expDim, 1]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [expDim, 1]);

% Check if matrix is symmetric.
verifyEqual(testCase, U, U', 'AbsTol', 1e-15);

end