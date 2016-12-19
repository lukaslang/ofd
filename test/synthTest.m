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
function tests = synthTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create triangulation of unit sphere.
[Faces, Verts] = sphTriang(3);
n = size(Faces, 1);

% Create vector spherical harmonics.
N = 1:10;
[Y1, Y2] = vspharmn(N, Faces, Verts);

% Create coefficients.
u = zeros(2*(N(end)^2 + 2*N(end)), 1);
% Compute vector spherical harmonics synthesis.
U = synth(cat(2, Y1, Y2), u);

% Check results.
verifyFalse(testCase, isempty(U));
verifyEqual(testCase, size(U), [n, 3]);
verifyEqual(testCase, U, zeros(n, 3));

end