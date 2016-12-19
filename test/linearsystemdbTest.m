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
function tests = linearsystemdbTest
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
n = size(Faces, 1);

% Create random data.
f1 = randi(255, m, 1);
f2 = randi(255, m, 1);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
M = 5;
N = 5;
[dim1, dim2, U, V, W, d1, d2, b] = linearsystemdb(Faces, Verts, 1:M, 1:N, f1, f2, h, tol);

% Check results.
verifyEqual(testCase, dim1 + dim2, 4*(N^2 + 2*N));
verifyFalse(testCase, isempty(U));
verifyEqual(testCase, size(U), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(V), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(W));
verifyEqual(testCase, size(W), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(d1));
verifyTrue(testCase, isvector(d1));
verifyEqual(testCase, size(d1), [2*(N^2 + 2*N), 1]);
verifyFalse(testCase, isempty(d2));
verifyTrue(testCase, isvector(d2));
verifyEqual(testCase, size(d2), [2*(N^2 + 2*N), 1]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [4*(N^2 + 2*N), 1]);

% Check if matrices are equal.
verifyEqual(testCase, U, V, 'AbsTol', 1e-15);
verifyEqual(testCase, U, W);

% Try with intervals.
M = 1:5;
N = 1:5;
[dim1i, dim2i, Ui, Vi, Wi, d1i, d2i, bi] = linearsystemdb(Faces, Verts, M, N, f1, f2, h, tol);
verifyEqual(testCase, dim1, dim1i, 'AbsTol', 1e-15);
verifyEqual(testCase, dim2, dim2i, 'AbsTol', 1e-15);
verifyEqual(testCase, Ui, U, 'AbsTol', 1e-15);
verifyEqual(testCase, Vi, V, 'AbsTol', 1e-15);
verifyEqual(testCase, Wi, W, 'AbsTol', 1e-15);
verifyEqual(testCase, d1i, d1, 'AbsTol', 1e-15);
verifyEqual(testCase, d2i, d2, 'AbsTol', 1e-15);
verifyEqual(testCase, bi, b, 'AbsTol', 1e-15);

end

function differentIntervalsTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(V, 1);
n = size(F, 1);

% Create random data.
f1 = randi(1, m, 1);
f2 = randi(1, m, 1);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
M = 1:3;
N = 4:5;
[dim1, dim2, U, V, W, d1, d2, b] = linearsystemdb(F, V, M, N, f1, f2, h, tol);

% Check results.
verifyEqual(testCase, dim1, 2*(M(end)^2 + 2*M(end) - M(1)^2 + 1));
verifyEqual(testCase, dim2, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1));
verifyFalse(testCase, isempty(U));
verifyEqual(testCase, size(U), [dim1, dim1]);
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(V), [dim2, dim2]);
verifyFalse(testCase, isempty(W));
verifyEqual(testCase, size(W), [dim1, dim2]);
verifyFalse(testCase, isempty(d1));
verifyTrue(testCase, isvector(d1));
verifyEqual(testCase, size(d1), [dim1, 1]);
verifyFalse(testCase, isempty(d2));
verifyTrue(testCase, isvector(d2));
verifyEqual(testCase, size(d2), [dim2, 1]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [dim1 + dim2, 1]);

% Check if matrices are symmetric.
verifyEqual(testCase, U', U, 'AbsTol', 1e-15);
verifyEqual(testCase, V', V, 'AbsTol', 1e-15);

end