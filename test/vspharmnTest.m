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
function tests = vspharmnTest
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

% Create vector spherical harmonics up to degree N.
N = 5;
[Y1, Y2] = vspharmn(1:N, F, V);
verifyFalse(testCase, isempty(Y1));
verifyFalse(testCase, isempty(Y2));
verifyEqual(testCase, size(Y1), [n, (N^2 + 2*N), 3]);
verifyEqual(testCase, size(Y2), [n, (N^2 + 2*N), 3]);

end

function invalidIntervalTest(testCase)
try
    vspharm([], [], []);
    verifyTrue(testCase, false, 'Assertion is required to fail.');  
end
try
    vspharm(0:3, [], []);
    verifyTrue(testCase, false, 'Assertion is required to fail.');  
end
try
    vspharm(-1:3, [], []);
    verifyTrue(testCase, false, 'Assertion is required to fail.');  
end   
try
    vspharm([1, 1, 2], [], []);
    verifyTrue(testCase, false, 'Assertion is required to fail.');  
end
try
    vspharm([1, 2, 4], [], []);
    verifyTrue(testCase, false, 'Assertion is required to fail.');  
end
try
    vspharm([1, 3, 2], [], []);
    verifyTrue(testCase, false, 'Assertion is required to fail.');  
end
try
    vspharm([3, 2, 1], [], []);
    verifyTrue(testCase, false, 'Assertion is required to fail.');  
end
end

function intervalTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

% Create vector spherical harmonics up to degree N.
l = 2;
Nu = 5;
[Y1, Y2] = vspharmn(l:Nu, F, V);
verifyFalse(testCase, isempty(Y1));
verifyFalse(testCase, isempty(Y2));
verifyEqual(testCase, size(Y1), [n, (Nu^2 + 2*Nu - l^2 + 1), 3]);
verifyEqual(testCase, size(Y2), [n, (Nu^2 + 2*Nu - l^2 + 1), 3]);

end

function interval2Test(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);

% Create vector spherical harmonics up to degree 5.
[Y1, Y2] = vspharmn(5, F, V);
verifyEqual(testCase, size(Y1, 2), 11);
verifyEqual(testCase, size(Y2, 2), 11);

[Y1, Y2] = vspharmn(1:5, F, V);
verifyEqual(testCase, size(Y1, 2), 35);
verifyEqual(testCase, size(Y2, 2), 35);

end

function onbTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);

% Create vector spherical harmonics up to degree N.
N = 1:5;
[Y1, Y2] = vspharmn(N, F, V);
Y = cat(2, Y1, Y2);

% Compute triangle areas.
a = triangArea(F, V);

% Check if surface integral over vector spherical harmonics are zero/one.
for k=1:size(Y, 2)
    for l=1:size(Y, 2)
        % Compute surface integral.
        v = triangIntegral(F, V, dot(Y(:, k, :), Y(:, l, :), 3), a);
        verifyTrue(testCase, isscalar(v));
        if(l == k)
            verifyEqual(testCase, v, 1, 'AbsTol', 1e-2);
        else
            verifyEqual(testCase, v, 0, 'AbsTol', 1e-2);
        end
    end
end
end