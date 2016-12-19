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
function tests = projecttoplaneTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

v = [pi, -pi, pi/2];
p = projecttoplane(v);
verifyFalse(testCase, isempty(p));
verifyEqual(testCase, p(:, 3), 0);
verifyEqual(testCase, sqrt(sum(p.^2, 2)), 3*pi/2, 'AbsTol', 1e-15);

end

function zeroTest(testCase)

v = [0, 0, 0];
p = projecttoplane(v);
verifyFalse(testCase, isempty(p));
verifyEqual(testCase, p, [0, 0, 0]);

end