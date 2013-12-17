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

% Define dataset and get result files.
name = 'cxcr4aMO2_290112';
resultsPath = fullfile('./', 'results', name, 'ofd');
resultsname = '2013-12-04-20-58-37-frames-114-116-unfiltered-1-100-7';
load(fullfile(resultsPath, sprintf('%s.mat', resultsname)));

% Import data.
disp('Loading precomputed data.');
path = fullfile('./', 'data', name, 'generated');
filename = 'frames-114-116-unfiltered-1-100-7';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));

% Load colormap for proper visualisation.
load(fullfile('./', 'data', name, 'cmapblue.mat'));

renderPath = fullfile('./', 'renderings', name, 'ofd', resultsname);
mkdir(renderPath);
mkdir(fullfile(renderPath, 'detail2'));

% Plot streamlines.
n = size(D.Faces, 1);
T = TriRep(D.Faces, D.Verts);
P = T.incenters;

% Select experiments.
idx = 208;

% Select section of cell division.
lim = [-0.2348, -0.1382, 0.1212, 0.2178];

% Vector length scaling.
scale = 0;
% Line width.
lw = 1;
% Face alpha.
fa = 0.25;

for k=idx

% Compute flow and scaling.
U = E{k}.U1 + E{k}.U2;
V = E{k}.V1 + E{k}.V2;
W = U + V;
nmax = max(sqrt(sum(W.^2, 2)));

G = figure;
axis square;
daspect([1, 1, 1]);
hold on;
colormap(cmap);
set(G, 'renderer', 'opengl');

trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), D.f{1}, 'EdgeColor', 'none', 'FaceAlpha', fa);
shading interp;
quiver3(P(:, 1), P(:, 2), P(:, 3), W(:, 1), W(:, 2), W(:, 3), scale, 'm', 'LineWidth', lw);
view(2);
axis(lim);
adjustFigure3;

H = figure;
axis square;
daspect([1, 1, 1]);
hold on;
colormap(cmap);
set(H, 'renderer', 'opengl');

trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), D.f{1}, 'EdgeColor', 'none', 'FaceAlpha', fa);
shading interp;
quiver3(P(:, 1), P(:, 2), P(:, 3), U(:, 1), U(:, 2), U(:, 3), scale, 'r', 'LineWidth', lw);
quiver3(P(:, 1), P(:, 2), P(:, 3), V(:, 1), V(:, 2), V(:, 3), scale, 'k', 'LineWidth', lw);
view(2);
axis(lim);
adjustFigure3;

I = figure;
axis square;
daspect([1, 1, 1]);
hold on;
colormap(cmap);
set(I, 'renderer', 'opengl');

trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), D.f{2}, 'EdgeColor', 'none', 'FaceAlpha', fa);
shading interp;
view(2);
axis(lim);
adjustFigure3;

% Save figures.
export_fig(fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-1.png', filename, k)), '-png', '-r300', '-transparent', G);
export_fig(fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-2.png', filename, k)), '-png', '-r300', '-transparent', H);
export_fig(fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-3.png', filename, k)), '-png', '-r300', '-transparent', I);

end