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
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
% Set working directory.
path = fullfile('./', 'data', name);

% Import data.
disp('Loading image data.');
%load(fullfile(path, 'frames-114-116-filtered.mat'));
%load(fullfile(path, 'frames-114-116-unfiltered.mat'));
load(fullfile(path, 'frames-140-142-unfiltered.mat'));

% Import cell centres.
disp('Loading cell centres.');
load(fullfile(path, 'thresholdedcenters.mat'));

% Load colormap for proper visualisation.
load(fullfile(path, 'cmapblue.mat'));

%frame = 114;
frame = 140;

% Set parameters.
Ns = 0:10;
alpha = 0.25;
s = 1;

% Prepare cell centres.
DX = F{frame}.X;
DY = F{frame}.Y;
DZ = -4.2832 * F{frame}.Z;
shift = -min(DZ);

% Fit sphere.
sc = mean([DX, DY, DZ + shift]);
sr = 300;
[sc, sr] = spherefit([DX, DY, DZ + shift], sc, sr);
DX = DX - sc(1);
DY = DY - sc(2);
DZ = DZ - sc(3) + shift;

% Fit spherical surface.
c = surffit(Ns, [DX, DY, DZ], alpha, s);

% Create triangulation.
[F, V] = sphTriang(7);

% Compute synthesis at vertices.
[Vs, rho] = surfsynth(Ns, V, c);

% Plot function rho on the unit sphere.
figure;
hold on;
axis([-1, 1, -1, 1, -1, 1]);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), rho, 'EdgeColor', 'none');
shading interp;
daspect([1, 1, 1]);
view(3);
colorbar;

% Plot surface.
figure;
hold on;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3));
shading interp;
daspect([1, 1, 1]);
view(3);

figure;
f = cell(2);
for k=1:2
    % Prepare data.
    img = flipdim(U{k}.u, 3);

    % Project data.
    [um, un, uo] = size(img);
    [X, Y, Z] = ndgrid(1:um, 1:un, 1:uo);

    % Compute radial maximum intensity projection.
    rs = linspace(0.8, 1.2, 80);
    VB = kron(rs', Vs);
    fb = dataFromCube(sc(1)+VB(:, 1), sc(2)+VB(:, 2), sc(3)+VB(:, 3), X, Y, 4.2832 * Z, img);
    f{k} = max(reshape(fb, size(Vs, 1), length(rs)), [], 2);

    subplot(2, 2, k);
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f{k}, 'EdgeColor', 'none');
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    
    subplot(2, 2, 2+k);
    axis([-1, 1, -1, 1, -1, 1]);
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f{k}, 'EdgeColor', 'none');
    shading interp;
    daspect([1, 1, 1]);
    view(3);
end

% Free memory.
clear VB;
clear U;
clear X;
clear Y;
clear Z;
clear fb;

% Set parameters.
N = 10;
h = 1;
alpha = 1;

disp('Computing optical flow...');
tic;
u = ofs(N, Ns, c, F, V, f{1}, f{2}, h, alpha);
toc;

% Get incenters of triangles.
TR = TriRep(F, Vs);
IC = TR.incenters;

% Plot result.
figure;
hold on;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f{1});
shading interp;
colormap(cmap);
daspect([1, 1, 1]);
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), u(:, 1), u(:, 2), u(:, 3), 0, 'r');

% Project and scale flow.
up = projecttoplane(u);

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
c = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

figure;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);