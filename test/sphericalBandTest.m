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
function test_suite = sphericalBandTest
    initTestSuite;
end

function resultTest

m = 10;
n = 10;
k = 5;
[x, y, z] = sphericalBand(m, n, linspace(5, 10, k));

assertFalse(isempty(x));
assertFalse(isempty(y));
assertFalse(isempty(z));

assertEqual(size(x), [m, n, k]);
assertEqual(size(y), [m, n, k]);
assertEqual(size(z), [m, n, k]);

end

function singleSphereTest

m = 10;
n = 10;
r = 1;
[x, y, z] = sphericalBand(m, n, r);

assertFalse(isempty(x));
assertFalse(isempty(y));
assertFalse(isempty(z));

assertEqual(size(x), [m, n]);
assertEqual(size(y), [m, n]);
assertEqual(size(z), [m, n]);

end