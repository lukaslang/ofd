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
function test_suite = vspharmnTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

% Create vector spherical harmonics up to degree N.
N = 5;
[Y, d] = vspharmn(N, F, V);
assertFalse(isempty(Y));
assertEqual(size(Y), [n, 2*(N^2 + 2*N), 3]);
assertFalse(isempty(d));
assertTrue(isvector(d));
assertEqual(size(d), [2*(N^2 + 2*N), 1]);

end

function invalidIntervalTest
try
    vspharm([], [], []);
    assertTrue(false, 'Assertion is required to fail.');  
end
try
    vspharm(0:3, [], []);
    assertTrue(false, 'Assertion is required to fail.');  
end
try
    vspharm(-1:3, [], []);
    assertTrue(false, 'Assertion is required to fail.');  
end   
try
    vspharm([1, 1, 2], [], []);
    assertTrue(false, 'Assertion is required to fail.');  
end
try
    vspharm([1, 2, 4], [], []);
    assertTrue(false, 'Assertion is required to fail.');  
end
try
    vspharm([1, 3, 2], [], []);
    assertTrue(false, 'Assertion is required to fail.');  
end
try
    vspharm([3, 2, 1], [], []);
    assertTrue(false, 'Assertion is required to fail.');  
end
end

function intervalTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

% Create vector spherical harmonics up to degree N.
l = 2;
Nu = 5;
[Y, d] = vspharmn(l:Nu, F, V);
assertFalse(isempty(Y));
assertEqual(size(Y), [n, 2*(Nu^2 + 2*Nu - l^2 + 1), 3]);
assertFalse(isempty(d));
assertTrue(isvector(d));
assertEqual(size(d), [2*(Nu^2 + 2*Nu - l^2 + 1), 1]);

end

function interval2Test

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);

% Create vector spherical harmonics up to degree 5.
[Y1, d1] = vspharmn(5, F, V);
[Y2, d2] = vspharmn(1:5, F, V);
assertAlmostEqual(Y1, Y2);
assertAlmostEqual(d1, d2);

end

function onbTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);

% Create vector spherical harmonics up to degree N.
N = 5;
[Y, ~] = vspharmn(N, F, V);

% Compute triangle areas.
a = triangArea(F, V);

% Check if surface integral over vector spherical harmonics are zero/one.
for k=1:size(Y, 2)
    for l=1:size(Y, 2)
        % Compute surface integral.
        v = triangIntegral(F, V, dot(Y(:, k, :), Y(:, l, :), 3), a);
        assertTrue(isscalar(v));
        if(l == k)
            assertAlmostEqual(v, 1, 1e-2);
        else
            assertAlmostEqual(v, 0, 1e-2);
        end
    end
end
end