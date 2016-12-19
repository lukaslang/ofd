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
function tests = residualTest
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
[F, V] = sphTriang(3);
m = size(V, 1);
n = size(F, 1);

% Create two images.
f1 = zeros(m, 1);
f2 = zeros(m, 1);

N = 5;
h = 1;
alpha = 1;

u = of(N, F, V, f1, f2, h, alpha);

verifyFalse(testCase, isempty(u));
verifyEqual(testCase, size(u), [n, 3]);
verifyEqual(testCase, u, zeros(n, 3));

[res, min] = residual(u, F, V, f1, f2, 1e-6);
verifyFalse(testCase, isempty(res));
verifyFalse(testCase, isempty(min));
verifyEqual(testCase, res, 0, 'AbsTol', 1e-15);
verifyEqual(testCase, min, 0, 'AbsTol', 1e-15);

end

function ofTest(testCase)

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

N = 3;
h = 1;
alpha = 1;

u = of(N, F, V, f1, f2, h, alpha);
verifyFalse(testCase, isempty(u));
verifyEqual(testCase, size(u), [n, 3]);

tic;
[res, min] = residual(u, F, V, f1, f2, 1e-6);
toc;
verifyTrue(testCase, res > 0);
verifyEqual(testCase, min, 0, 'AbsTol', 1e-15);

end