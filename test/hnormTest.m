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
function test_suite = hnormTest
    initTestSuite;
end

function resultTest

s = 1;
N = 1:10;
a = vspharmeigs(N);
c = zeros(length(a), 1);

% Compute Sobolev norm.
n = hnorm(s, a, c);
assertTrue(isscalar(n));
assertFalse(isempty(n));
assertAlmostEqual(n, 0);

end

function ofTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Create two images.
Ynj = spharm(5, V);
f1 = Ynj(:, 3);
% Create rotation matrix.
theta = -pi/9;
T = [   cos(theta),     sin(theta), 0;
       -sin(theta),     cos(theta), 0;
                 0,              0, 1];
% Create rotated image.
Ynj = spharm(5, V*T);
f2 = Ynj(:, 3);

N = 5;
s = -1;
h = 1;
alpha = 10;

% Compute functions needed for solving the linear system.
[dim, U, d, b] = linearsystem(F, V, 1:N, f1, f2, h, 1e-6);

% Solve linear system.
[u, ~] = ofsolve(dim, U, b, d, alpha, s, 30);
assertFalse(isempty(u));
assertTrue(isvector(u));
assertEqual(length(u), 2*(N^2+2*N));

% Compute Sobolev norm of the oscillating result.
n = hnorm(s, d, u);
assertTrue(isscalar(n));
assertFalse(isempty(n));
assertTrue(n > 0);

% Compute Sobolev norm of space H^1. The recovered solution should have a
% larger norm in the space H^1.
ns = hnorm(-s, d, u);
assertTrue(ns > n);

end