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
function test_suite = ofdbsolveTest
    initTestSuite;
end

function sameBasisTest

% Create triangulation of unit sphere.
[Faces, Verts] = sphTriang(3);
m = size(Verts, 1);
n = size(Faces, 1);

% Create two images.
f1 = randi(255, m, 1);
f2 = randi(255, m, 1);

M = 5;
N = 5;
h = 1;
alpha = 1;
beta = 1;

% Compute linear system.
[dim1, dim2, U, V, W, d1, d2, Y1, Y2, bi] = linearsystemdb(Faces, Verts, M, N, f1, f2, h, 1e-6);
% Solve linear system.
[Ui, Vi, ui, vi, ~] = ofdbsolve(dim1, dim2, U, V, W, d1, d2, Y1, Y2, bi, alpha, beta, 1, -1);

% Compute functions for same basis.
[dim, At, d, Y, b] = linearsystem(Faces, Verts, N, f1, f2, h, 1e-6);
assertAlmostEqual(dim, dim1);
assertAlmostEqual(dim, dim2);
assertAlmostEqual(At, U);
assertAlmostEqual(At, V);
assertAlmostEqual(At, W);
assertAlmostEqual(d, d1);
assertAlmostEqual(d, d2);
assertAlmostEqual(Y, Y1);
assertAlmostEqual(Y, Y2);
assertAlmostEqual([b; b], bi);

% Solve linear system.
[U, V, u, v, ~] = ofdsolve(dim, At, b, Y, d, alpha, beta, 1, -1);
assertAlmostEqual(U, Ui, 1e-10);
assertAlmostEqual(V, Vi, 1e-10);
assertAlmostEqual(u, ui, 1e-10);
assertAlmostEqual(v, vi, 1e-10);

end