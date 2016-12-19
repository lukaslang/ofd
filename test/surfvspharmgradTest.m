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
function tests = surfvspharmgradTest
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
[F, V] = sphTriang(4);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Compute orthonormal basis of the tangent space at incenters.
B = surftangentialbasis(Ns, c, F, V);
[~, E1, E2] = orthonormalise(squeeze(B(:, 1, :)), squeeze(B(:, 2, :)));

% Compute covariant derivative of vector spherical harmonics of degrees N.
N = 10;
[Z1, Z2, Z3, Z4] = surfvspharmgrad(E1, E2, F, V, N, Ns, c);
verifyFalse(testCase, isempty(Z1));
verifyFalse(testCase, isempty(Z2));
verifyFalse(testCase, isempty(Z3));
verifyFalse(testCase, isempty(Z4));
verifyEqual(testCase, size(Z1), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);
verifyEqual(testCase, size(Z2), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);
verifyEqual(testCase, size(Z3), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);
verifyEqual(testCase, size(Z4), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);

end