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
function tests = dataFromCubeTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

m = 30;
n = 20;
o = 10;

u = ones(m, n, o);

[x, y, z] = meshgrid(1:m, 1:n, 1:o);
[X, Y, Z] = ndgrid(1:m, 1:n, 1:o);
f = dataFromCube(x(:), y(:), z(:), X, Y, Z, u);

verifyFalse(testCase, isempty(f));
verifyEqual(testCase, size(f), [m*n*o, 1]);
verifyEqual(testCase, f, ones(m*n*o, 1));

end

function pointsOutsideCubeTest(testCase)

m = 30;
n = 20;
o = 10;

u = ones(m, n, o);

% Create data outside cube.
x = [m+1, 1, 1];
y = [1, n+1, 1];
z = [1, 1, o+1];

[X, Y, Z] = ndgrid(1:m, 1:n, 1:o);
f = dataFromCube(x(:), y(:), z(:), X, Y, Z, u);

verifyFalse(testCase, isempty(f));
verifyEqual(testCase, size(f), [3, 1]);
verifyEqual(testCase, f, zeros(3, 1));

end


function interpolatingTest(testCase)

m = 30;
n = 20;
o = 10;

u = ones(m, n, o);

[x, y, z] = meshgrid(1:m-1, 1:n-1, 1:o-1);

% Add fractional values.
x = x + 0.5;
y = y + 0.5;
z = z + 0.5;

[X, Y, Z] = ndgrid(1:m, 1:n, 1:o);
f = dataFromCube(x(:), y(:), z(:), X, Y, Z, u);

verifyFalse(testCase, isempty(f));
verifyEqual(testCase, size(f), [(m-1)*(n-1)*(o-1), 1]);
verifyEqual(testCase, f, ones((m-1)*(n-1)*(o-1), 1));

end