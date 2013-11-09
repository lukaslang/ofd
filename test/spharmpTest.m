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
function test_suite = spharmpTest
    initTestSuite;
end

function resultTest

% Create parametrisation and coordinates.
m = 10;
n = 20;

% Create parametrisation and polar coordinates.
[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));

% Create spherical harmonics.
Ynj = spharmp(5, phi, t);
assertFalse(isempty(Ynj));
assertEqual(size(Ynj), [m*n, 11]);

end

function visualiseTest

% Create parametrisation and coordinates.
m = 200;
n = 200;

% Create parametrisation and polar coordinates.
[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));

% Create spherical harmonics.
N = 3;
Ynj = spharmp(N, phi, t);
assertFalse(isempty(Ynj));
assertEqual(size(Ynj), [m*n, 2*N + 1]);

% Convert to cartesian coordinates.
x = sqrt(1 - t .^2) .* cos(phi);
y = sqrt(1 - t .^2) .* sin(phi);
z = t;

figure;
for k=1:2*N+1
    f = Ynj(:, k);
    f = reshape(f, m, n);

    % Visualise data.
    subplot(1, 2*N+1, k);
    surf(x, y, z, f, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
end

end