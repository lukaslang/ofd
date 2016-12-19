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
function tests = maximumIntensityTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create cube.
u = ones(100, 100, 100);

% Specify sphere in the interior of u.
r = 20;
c = [50, 50, 50];

% Specify parametrisation of band.
m = 20;
n = 20;
rs = linspace(r-5, r+5, 10);

[X, Y, Z] = ndgrid(1:100, 1:100, 1:100);
f = maximumIntensity(c, m, n, rs, X, Y, Z, u);
verifyFalse(testCase, isempty(f));
verifyEqual(testCase, f, ones(m, n));

end

function halfSphereTest(testCase)

% Create cube.
u = ones(100, 100, 100);

% Specify sphere with only upper half in u.
r = 20;
c = [50, 50, 0];

% Specify parametrisation of band.
m = 20;
n = 20;
rs = linspace(r-5, r+5, 10);

[X, Y, Z] = ndgrid(1:100, 1:100, 1:100);
f = maximumIntensity(c, m, n, rs, X, Y, Z, u);
verifyFalse(testCase, isempty(f));
verifyEqual(testCase, f, [zeros(m, n/2), ones(m, n/2)]);

end

function regularSphereTest(testCase)

% Create cube.
u = ones(100, 100, 100);

% Specify sphere in the interior of u.
c = [50, 50, 50];

% Specify parametrisation of single sphere.
m = 20;
n = 20;
rs = 20;

[X, Y, Z] = ndgrid(1:100, 1:100, 1:100);
f = maximumIntensity(c, m, n, rs, X, Y, Z, u);
verifyFalse(testCase, isempty(f));
verifyEqual(testCase, f, ones(m, n));

end