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
function test_suite = linearsystemdbTest
    initTestSuite;
end

function resultTest

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
assertEqual(dim1 + dim2, 4*(N^2 + 2*N));
assertFalse(isempty(U));
assertEqual(size(U), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
assertFalse(isempty(V));
assertEqual(size(V), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
assertFalse(isempty(W));
assertEqual(size(W), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
assertFalse(isempty(d1));
assertTrue(isvector(d1));
assertEqual(size(d1), [2*(N^2 + 2*N), 1]);
assertFalse(isempty(d2));
assertTrue(isvector(d2));
assertEqual(size(d2), [2*(N^2 + 2*N), 1]);
assertFalse(isempty(b));
assertTrue(isvector(b));
assertEqual(size(b), [4*(N^2 + 2*N), 1]);

% Check if matrices are equal.
assertAlmostEqual(U, V);
assertAlmostEqual(U, W);

% Try with intervals.
M = 1:5;
N = 1:5;
[dim1i, dim2i, Ui, Vi, Wi, d1i, d2i, bi] = linearsystemdb(Faces, Verts, M, N, f1, f2, h, tol);
assertAlmostEqual(dim1, dim1i);
assertAlmostEqual(dim2, dim2i);
assertAlmostEqual(Ui, U);
assertAlmostEqual(Vi, V);
assertAlmostEqual(Wi, W);
assertAlmostEqual(d1i, d1);
assertAlmostEqual(d2i, d2);
assertAlmostEqual(bi, b);

end

function differentIntervalsTest

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
assertEqual(dim1, 2*(M(end)^2 + 2*M(end) - M(1)^2 + 1));
assertEqual(dim2, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1));
assertFalse(isempty(U));
assertEqual(size(U), [dim1, dim1]);
assertFalse(isempty(V));
assertEqual(size(V), [dim2, dim2]);
assertFalse(isempty(W));
assertEqual(size(W), [dim1, dim2]);
assertFalse(isempty(d1));
assertTrue(isvector(d1));
assertEqual(size(d1), [dim1, 1]);
assertFalse(isempty(d2));
assertTrue(isvector(d2));
assertEqual(size(d2), [dim2, 1]);
assertFalse(isempty(b));
assertTrue(isvector(b));
assertEqual(size(b), [dim1 + dim2, 1]);

% Check if matrices are symmetric.
assertAlmostEqual(U', U);
assertAlmostEqual(V', V);

end