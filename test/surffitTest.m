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
function test_suite = surffitTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[~, V] = sphTriang(5);
n = size(V, 1);

% Set parameters.
N = 0:5;
s = 1;
alpha = 1;

% Fit surface to sphere.
[c, Y] = surffit(N, V, alpha, s);
assertFalse(isempty(c));
assertFalse(isempty(Y));

% Recover function at specified vertices.
rho = Y*c;

% Compute coordinates of fitted surface at data points.
F = bsxfun(@times, V, rho);

% Check if points are still on unit sphere.
len = sqrt(sum(F .^ 2,2));
assertAlmostEqual(len, ones(n, 1), 1e-12);

% Check if points equal.
assertAlmostEqual(F, V, 1e-12);

end

function result2Test

% Create triangulation of unit sphere.
[F, V] = sphTriang(6);

% Create data using one spherical harmonic of degree n and order m.
n = 5;
m = 3;
Y = spharm(n, V);
X = bsxfun(@times, V, 1 + Y(:, m));

% Plot surface.
figure;
hold on;
trisurf(F, X(:, 1), X(:, 2), X(:, 3));
shading interp;
daspect([1, 1, 1]);
view(3);

% Set parameters.
N = 0:5;
s = 1;
alpha = 1;

% Fit surface.
[c, Y] = surffit(N, X, alpha, s);
assertFalse(isempty(c));
assertFalse(isempty(Y));

% Recover function at specified vertices.
rho = Y*c;

% Compute coordinates of fitted surface at data points.
S = bsxfun(@times, V, rho);

% Check if radial error of fitted points is small.
lenS = sqrt(sum(S .^ 2,2));
lenX = sqrt(sum(X .^ 2,2));
assertAlmostEqual(lenS, lenX, 1e-2);

% Check if points are almost equal.
assertAlmostEqual(S, X, 1e-2);

% Plot surface.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3));
shading interp;
daspect([1, 1, 1]);
view(3);

end