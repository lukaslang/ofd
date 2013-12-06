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
function test_suite = projecttoplaneTest
    initTestSuite;
end

function resultTest

v = [pi, -pi, pi/2];
p = projecttoplane(v);
assertFalse(isempty(p));
assertEqual(p(:, 3), 0);
assertAlmostEqual(sqrt(sum(p.^2, 2)), 3*pi/2);

end