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
function test_suite = spharmTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[~, V] = sphTriang(3);

% Create spherical harmonics.
N = 5;
Ynj = spharm(N, V);
assertFalse(isempty(Ynj));
assertEqual(size(Ynj), [size(V, 1), 2*N + 1]);

end

function visualiseTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Create spherical harmonics.
N = 3;
Ynj = spharm(N, V);
assertFalse(isempty(Ynj));
assertEqual(size(Ynj), [size(V, 1), 2*N + 1]);

figure;
for k=1:2*N + 1
    % Visualise data.
    subplot(2, 2*N + 1, k);
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), Ynj(:, k), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    % Interpolated shading.
    subplot(2, 2*N + 1, k + 2*N + 1);
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), Ynj(:, k), 'EdgeColor', 'none');
    shading interp;
    daspect([1, 1, 1]);
    view(3);
end

end