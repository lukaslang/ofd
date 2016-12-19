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
function tests = vspharmsynthTest
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
n = size(Faces, 1);

N = 10;
% Create coefficients.
u = zeros(2*(N^2 + 2*N), 1);
% Compute vector spherical harmonics synthesis.
[U1, U2] = vspharmsynth(1:N, Faces, Verts, u);

% Check results.
verifyFalse(testCase, isempty(U1));
verifyFalse(testCase, isempty(U2));
verifyEqual(testCase, size(U1), [n, 3]);
verifyEqual(testCase, size(U2), [n, 3]);
verifyEqual(testCase, U1, zeros(n, 3));
verifyEqual(testCase, U2, zeros(n, 3));

end

function intervalTest(testCase)

% Create triangulation of unit sphere.
[Faces, Verts] = sphTriang(3);
n = size(Faces, 1);

N = 3:10;
% Create coefficients.
u = zeros(2*(N(end)^2 + 2*N(end) - N(1)^2 + 1), 1);
% Compute vector spherical harmonics synthesis.
[U1, U2] = vspharmsynth(N, Faces, Verts, u);

% Check results.
verifyFalse(testCase, isempty(U1));
verifyFalse(testCase, isempty(U2));
verifyEqual(testCase, size(U1), [n, 3]);
verifyEqual(testCase, size(U2), [n, 3]);
verifyEqual(testCase, U1, zeros(n, 3));
verifyEqual(testCase, U2, zeros(n, 3));

N = 1:10;
% Create coefficients.
u = zeros(2*(N(end)^2 + 2*N(end) - N(1)^2 + 1), 1);
% Compute vector spherical harmonics synthesis.
[U1, U2] = vspharmsynth(N, Faces, Verts, u);

% Check results.
verifyFalse(testCase, isempty(U1));
verifyFalse(testCase, isempty(U2));
verifyEqual(testCase, size(U1), [n, 3]);
verifyEqual(testCase, size(U2), [n, 3]);
verifyEqual(testCase, U1, zeros(n, 3));
verifyEqual(testCase, U2, zeros(n, 3));

end

function memConstraintTest(testCase)

% Create triangulation of unit sphere.
[Faces, Verts] = sphTriang(3);
n = size(Faces, 1);

N = 10;
% Create coefficients.
u = zeros(2*(N^2 + 2*N), 1);
% Compute vector spherical harmonics synthesis.
mem = 2e6;
[U1, U2] = vspharmsynth(1:N, Faces, Verts, u, mem);

% Check results.
verifyFalse(testCase, isempty(U1));
verifyFalse(testCase, isempty(U2));
verifyEqual(testCase, size(U1), [n, 3]);
verifyEqual(testCase, size(U2), [n, 3]);
verifyEqual(testCase, U1, zeros(n, 3));
verifyEqual(testCase, U2, zeros(n, 3));

end

function matrixTest(testCase)

% Create triangulation of unit sphere.
[Faces, Verts] = sphTriang(3);
n = size(Faces, 1);

N = 10;
% Create coefficient matrix.
dim = 5;
u = zeros(2*(N^2 + 2*N), dim);
% Compute vector spherical harmonics synthesis.
mem = 2e9;
[U1, U2] = vspharmsynth(1:N, Faces, Verts, u, mem);

% Check results.
verifyFalse(testCase, isempty(U1));
verifyFalse(testCase, isempty(U2));
verifyEqual(testCase, size(U1), [n, 3, dim]);
verifyEqual(testCase, size(U2), [n, 3, dim]);
verifyEqual(testCase, U1, zeros(n, 3, dim));
verifyEqual(testCase, U2, zeros(n, 3, dim));

end

function matrixPerformanceTest(testCase)

% Create triangulation of unit sphere.
[Faces, Verts] = sphTriang(3);

N = 30;
% Create coefficient matrix.
dim = 2;
u = zeros(2*(N^2 + 2*N), dim);
% Compute vector spherical harmonics synthesis.
mem = 2e9;
t = tic;
vspharmsynth(1:N, Faces, Verts, u, mem);
elapsed = toc(t);
disp(elapsed);

u = zeros(2*(N^2 + 2*N), 1);
t = tic;
for k=1:dim
    vspharmsynth(1:N, Faces, Verts, u, mem);
end
elapsed2 = toc(t);
disp(elapsed2);

verifyTrue(testCase, elapsed < elapsed2);

end