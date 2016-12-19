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
function tests = sphericalBandTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

m = 10;
n = 10;
k = 5;
[x, y, z] = sphericalBand(m, n, linspace(5, 10, k));

verifyFalse(testCase, isempty(x));
verifyFalse(testCase, isempty(y));
verifyFalse(testCase, isempty(z));

verifyEqual(testCase, size(x), [m, n, k]);
verifyEqual(testCase, size(y), [m, n, k]);
verifyEqual(testCase, size(z), [m, n, k]);

end

function singleSphereTest(testCase)

m = 10;
n = 10;
r = 1;
[x, y, z] = sphericalBand(m, n, r);

verifyFalse(testCase, isempty(x));
verifyFalse(testCase, isempty(y));
verifyFalse(testCase, isempty(z));

verifyEqual(testCase, size(x), [m, n]);
verifyEqual(testCase, size(y), [m, n]);
verifyEqual(testCase, size(z), [m, n]);

end