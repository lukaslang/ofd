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
function test_suite = orthonormaliseTest
    initTestSuite;
end

function resultTest

% Create two test vectors.
v = [3, 1, 0];
w = [2, 2, 0];

% Orthonormalise.
[A, e, f] = orthonormalise(v, w);

% Check first basis vector.
assertAlmostEqual(e, [3/sqrt(10), 1/sqrt(10), 0]);
% Check second basis vector.
assertAlmostEqual(f, [-1/sqrt(10), 3/sqrt(10), 0]);

% Check lengths.
assertAlmostEqual(sqrt(sum(e .^ 2, 2)), 1);
assertAlmostEqual(sqrt(sum(f .^ 2, 2)), 1);

% Check orthogonality.
assertAlmostEqual(e * f', 0);

% Check coefficients.
assertAlmostEqual(A(1, 1, :), 1/sqrt(10));
assertAlmostEqual(A(1, 2, :), 0);
assertAlmostEqual(A(2, 1, :), -8/(10 * sqrt(4/25 + 36/25)));
assertAlmostEqual(A(2, 2, :), 1/sqrt(4/25 + 36/25));

% Check that A transforms {v, w} into {e, f}.
assertAlmostEqual(A * [v; w], [e; f]);

end