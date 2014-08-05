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
function test_suite = ofsTest
    initTestSuite;
end

function identityMapTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
m = size(V, 1);
n = size(F, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Create two images.
f1 = zeros(m, 1);
f2 = zeros(m, 1);

N = 10;
h = 1;
alpha = 1;

u = ofs(N, Ns, c, F, V, f1, f2, h, alpha);

assertFalse(isempty(u));
assertEqual(size(u), [n, 3]);
assertEqual(u, zeros(n, 3));

end

function visualiseTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Create two images.
Ynj = spharm(5, V);
f1 = Ynj(:, 3);
% Create rotation matrix.
theta = -pi/9;
T = [   cos(theta),     sin(theta), 0;
       -sin(theta),     cos(theta), 0;
                 0,              0, 1];
% Create rotated image.
Ynj = spharm(5, V*T);
f2 = Ynj(:, 3);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

N = 5;
h = 1;
alpha = 1;

u = ofs(N, Ns, c, F, V, f1, f2, h, alpha);
assertFalse(isempty(u));
assertEqual(size(u), [n, 3]);

% Compute residual.
[res, ~] = residual(u, F, V, f1, f2, 1e-6);
fprintf('Residual: %f.\n', res);

% Compare to other optical flow implementation.
v = of(N, F, V, f1, f2, h, alpha);
assertFalse(isempty(v));
assertEqual(size(v), [n, 3]);
% Compute residual.
[res, ~] = residual(v, F, V, f1, f2, 1e-6);
fprintf('Residual: %f.\n', res);

% Check if flow is almost equal.
assertAlmostEqual(u, v, 0.05);

TR = TriRep(F, V);
P = TR.incenters;

figure;
subplot(1, 3, 1);
axis([-1, 1, -1, 1, -1, 1]);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), f1);
shading interp;
daspect([1, 1, 1]);
view(3);
subplot(1, 3, 2);
axis([-1, 1, -1, 1, -1, 1]);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), f2);
shading interp;
daspect([1, 1, 1]);
view(3);
subplot(1, 3, 3);
axis([-1, 1, -1, 1, -1, 1]);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), f1);
shading interp;
daspect([1, 1, 1]);
view(3);
quiver3(P(:, 1), P(:, 2), P(:, 3), u(:, 1), u(:, 2), u(:, 3), 0, 'm');
quiver3(P(:, 1), P(:, 2), P(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'b');

end