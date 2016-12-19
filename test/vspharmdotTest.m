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
function tests = vspharmdotTest
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
n = size(F, 1);

% Compute triangle areas.
a = triangArea(F, V);

% Arbitrary vector spherical harmonics.
deg = 1;
ord = 1;
[Y1, Y2] = vspharm(deg, F, V);
Y1nj = squeeze(Y1(:, ord, :));
Y2nj = squeeze(Y2(:, ord, :));

% Compute dot product with vector spherical harmonics of degrees N.
N = 1;
Z = vspharmdot(Y1nj, F, V, N);
verifyFalse(testCase, isempty(Z));
verifyEqual(testCase, size(Z), [n, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);
verifyEqual(testCase, triangIntegral(F, V, Z(:, 1), a), 1, 'AbsTol', 1e-2);
for k=2:size(Z, 2)
    verifyEqual(testCase, triangIntegral(F, V, Z(:, k), a), 0, 'AbsTol', 1e-2);
end

% Compute dot product with vector spherical harmonics of degrees N.
Z = vspharmdot(Y2nj, F, V, N);
verifyFalse(testCase, isempty(Z));
verifyEqual(testCase, size(Z), [n, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);
verifyEqual(testCase, triangIntegral(F, V, Z(:, 4), a), 1, 'AbsTol', 1e-2);
for k=[1:3, 5:6]
    verifyEqual(testCase, triangIntegral(F, V, Z(:, k), a), 0, 'AbsTol', 1e-2);
end

end