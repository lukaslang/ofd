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
function tests = spharmnTest
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
[~, V] = sphTriang(3);
m = size(V, 1);

% Create scalar spherical harmonics up to degree N.
N = 5;
Y = spharmn(0:N, V);
verifyFalse(testCase, isempty(Y));
verifyEqual(testCase, size(Y), [m, N^2 + 2*N + 1]);

end

function intervalTest(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(3);
m = size(V, 1);

% Create scalar spherical harmonics up to degree N.
l = 2;
Nu = 5;
Y = spharmn(l:Nu, V);
verifyFalse(testCase, isempty(Y));
verifyEqual(testCase, size(Y), [m, (Nu^2 + 2*Nu - l^2 + 1)]);

end

function interval2Test(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(3);

% Create scalar spherical harmonics up to degree 5.
Y = spharmn(5, V);
verifyEqual(testCase, size(Y, 2), 11);

Y = spharmn(1:5, V);
verifyEqual(testCase, size(Y, 2), 35);

end