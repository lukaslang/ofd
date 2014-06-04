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
function test_suite = surfsynthTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[~, V] = sphTriang(5);
n = size(V, 1);

% Set parameters.
N = 0;
% Compute coefficient such that rho is identically one, i.e. unit sphere.
c = 1/spharm(0, [0, 1, 0]);

% Compute synthesis.
[S, rho] = surfsynth(N, V, c);

% Check if points are still on unit sphere.
len = sqrt(sum(S .^ 2, 2));
assertAlmostEqual(len, ones(n, 1));

% Check if points equal.
assertAlmostEqual(S, V);

% Check if rho is identically one.
assertAlmostEqual(rho, ones(n, 1));

end

function resultIntervalTest

% Create triangulation of unit sphere.
[~, V] = sphTriang(5);
n = size(V, 1);

% Set parameters.
N = 0:5;
dim = N(end)^2 + 2*N(end) - N(1)^2 + 1;
% Compute coefficient such that rho is identically one, i.e. unit sphere.
c = [1/spharm(0, [0, 1, 0]); zeros(dim - 1, 1)];

% Compute synthesis.
[S, rho] = surfsynth(N, V, c);

% Check if points are still on unit sphere.
len = sqrt(sum(S .^ 2, 2));
assertAlmostEqual(len, ones(n, 1));

% Check if points equal.
assertAlmostEqual(S, V);

% Check if rho is identically one.
assertAlmostEqual(rho, ones(n, 1));

end

function sphereTest

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

% Compute synthesis from scratch.
[S, rho1] = surfsynth(N, V, c);

% Compute synthesis with Y.
[SY, rho2] = surfsynth(N, V, c, Y);
assertAlmostEqual(S, SY);
assertAlmostEqual(rho1, rho2);

% Check if points are still on unit sphere.
len = sqrt(sum(S .^ 2, 2));
assertAlmostEqual(len, ones(n, 1), 1e-12);

% Check if points equal.
assertAlmostEqual(S, V, 1e-12);

% Check if functions rho are identically one.
assertAlmostEqual(rho1, ones(n, 1), 1e-12);
assertAlmostEqual(rho2, ones(n, 1), 1e-12);

end