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
load(fullfile(path, 'frames-114-116-unfiltered.mat'));

% Import cell centres.
disp('Loading cell centres.');
load(fullfile(path, 'thresholdedcenters.mat'));

% Load colormap for proper visualisation.
load(fullfile(path, 'cmapblue.mat'));

frame = 114;

% Set decomposition parameters.
M = 1:5;
N = 6:10;
h = 1;
alpha = 1;
beta = 100;

% Prepare cell centres.
X = F{frame}.X;
Y = F{frame}.Y;
Z = -4.2832 * F{frame}.Z;
shift = -min(Z);

% Fit sphere.
sc = mean([X(:), Y(:), Z(:) + shift]);
sr = 300;
[c, r] = spherefit([X(:), Y(:), Z(:) + shift], sc, sr);

% Create triangulation of northern hemisphere of the fitted sphere.
[F, V] = halfsphTriang(7);

figure;
f = cell(2);
for k=1:2
    % Prepare data.
    u = flipdim(U{k}.u, 3);

    % Project data.
    [um, un, uo] = size(u);
    [X, Y, Z] = ndgrid(1:um, 1:un, 1:uo);

    % Compute radial maximum intensity projection.
    rs = linspace(r-40, r+40, 80);
    VB = kron(rs', V);
    fb = dataFromCube(c(1)+VB(:, 1), c(2)+VB(:, 2), c(3)+VB(:, 3), X, Y, 4.2832 * Z, u);
    f{k} = max(reshape(fb, size(V, 1), length(rs)), [], 2);

    subplot(1, 2, k);
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

% Compute decomposition.
disp('Computing decomposition...');
tic;
[u, v] = ofdb(M, N, F, V, f{1}, f{2}, h, alpha, beta);
toc;

% Get incenters of triangles.
TR = TriRep(F, V);
P = TR.incenters;

% Compute colour space scaling.
nmax = max(sqrt(sum((u + v).^2, 2)));

% Plot colour images of projections.
c = double(squeeze(computeColour(u(:, 1)/nmax, u(:, 2)/nmax))) ./ 255;
figure;
axis([-1, 1, -1, 1, 0, 1]);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

c = double(squeeze(computeColour(v(:, 1)/nmax, v(:, 2)/nmax))) ./ 255;
figure;
axis([-1, 1, -1, 1, 0, 1]);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

c = double(squeeze(computeColour((u(:, 1) + v(:, 1))/nmax, (u(:, 2) + v(:, 2))/nmax))) ./ 255;
figure;
axis([-1, 1, -1, 1, 0, 1]);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);