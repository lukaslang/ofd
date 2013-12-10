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
function test_suite = linearsystemTest
    initTestSuite;
end

function resultTest

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
assertEqual(dim, 2*(N^2 + 2*N));
assertFalse(isempty(U));
assertEqual(size(U), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
assertFalse(isempty(d));
assertTrue(isvector(d));
assertEqual(size(d), [2*(N^2 + 2*N), 1]);
assertFalse(isempty(b));
assertTrue(isvector(b));
assertEqual(size(b), [2*(N^2 + 2*N), 1]);

% Check if matrix is symmetric.
assertAlmostEqual(U, U');

end

function intervalTest

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
assertEqual(dim, expDim);
assertFalse(isempty(U));
assertEqual(size(U), [expDim, expDim]);
assertFalse(isempty(d));
assertTrue(isvector(d));
assertEqual(size(d), [expDim, 1]);
assertFalse(isempty(b));
assertTrue(isvector(b));
assertEqual(size(b), [expDim, 1]);

% Check if matrix is symmetric.
assertAlmostEqual(U, U');

end