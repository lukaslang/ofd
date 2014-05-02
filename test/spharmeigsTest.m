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
function test_suite = spharmeigsTest
    initTestSuite;
end

function resultTest

N = 0:3;
dim = N(end)^2 + 2*N(end) - N(1)^2 + 1;

% Compute eigenvalues.
d = spharmeigs(N);
assertTrue(isvector(d));
assertEqual(length(d), dim);

assertEqual(d(1), 0);
assertEqual(d(2:4), 2*ones(3, 1));
assertEqual(d(5:9), 6*ones(5, 1));
assertEqual(d(10:16), 12*ones(7, 1));

end

function intervalTest

N = 2:3;
dim = N(end)^2 + 2*N(end) - N(1)^2 + 1;

% Compute eigenvalues.
d = spharmeigs(N);
assertTrue(isvector(d));
assertEqual(length(d), dim);

assertEqual(d(1:5), 6*ones(5, 1));
assertEqual(d(6:12), 12*ones(7, 1));

end