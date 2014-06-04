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
N = 0:10;
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
c = surffit(N, [DX, DY, DZ], alpha, s);

% Create triangulation.
[F, V] = sphTriang(7);

% Compute synthesis at vertices.
[S, rho] = surfsynth(N, V, c);

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
trisurf(F, S(:, 1), S(:, 2), S(:, 3));
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
    VB = kron(rs', S);
    fb = dataFromCube(sc(1)+VB(:, 1), sc(2)+VB(:, 2), sc(3)+VB(:, 3), X, Y, 4.2832 * Z, img);
    f{k} = max(reshape(fb, size(S, 1), length(rs)), [], 2);

    subplot(2, 2, k);
    trisurf(F, S(:, 1), S(:, 2), S(:, 3), f{k}, 'EdgeColor', 'none');
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