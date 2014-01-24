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
function test_suite = triangsectionTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = halfsphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [8, 3]);
assertEqual(size(V), [8, 3]);
assert(all(V(:, 3) >= 0));

[Fs, Vs] = triangsection(F, V, [-1, 1, -1, 1, 0, 1]);

assertEqual(Fs, F);
assertEqual(Vs, V);

end

function visualiseTest

% Generate icosahedron.
tic;
[F, V] = halfsphTriang(4);
toc;
assertFalse(isempty(F));
assertFalse(isempty(V));
assert(all(V(:, 3) >= 0));

[Fs, Vs] = triangsection(F, V, [-1, 0, 0, 1, 0, 1]);

figure;
trisurf(Fs, Vs(:, 1), Vs(:, 2), Vs(:, 3));
daspect([1, 1, 1]);

end