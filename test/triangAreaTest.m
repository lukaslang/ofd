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
function test_suite = triangAreaTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Compute area of triangles.
A = triangArea(F, V);
assertFalse(isempty(A));
assertAlmostEqual(A, (4*pi/size(F, 1)) * ones(size(F, 1), 1), 1e-3);

end