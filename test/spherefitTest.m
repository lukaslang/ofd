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
function tests = spherefitTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang(3);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));

m = size(V, 1);

% Add noise to data.
V = V + 0.1 * randn(m, 3);

% Fit sphere.
[c, r] = spherefit(V, [0, 0, 0], 1);
fprintf('Fitted sphere is c=[%.4f, %.4f, %.4f], r=%.4f\n', c, r);

% Visualise fitted sphere.
figure;
hold on;
[X, Y, Z] = sphere(20);
surf(c(1) + r*X, c(2) + r*Y, c(3) + r*Z, 'FaceColor', 'b', 'FaceAlpha', 0.4);
daspect([1, 1, 1]);
scatter3(V(:, 1), V(:, 2), V(:, 3), 'r*');
view(3);

end