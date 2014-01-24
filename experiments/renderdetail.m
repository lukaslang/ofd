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

% Select experiments.
idx = 208;

% Select section of cell division.
lim = [-0.2348, -0.1382, 0.1212, 0.2178];

% Vector length scaling.
scale = 1;
% Line width.
lw = 1;
% Face alpha.
fa = 0.25;

% Restrict triangulation.
trlim = [-0.25, -0.13, 0.1, 0.23, 0, 1];
[Faces, Verts, Fidx, Vidx] = triangsection(D.Faces, D.Verts, trlim);
n = size(Faces, 1);
T = TriRep(Faces, Verts);
P = T.incenters;

for k=idx

% Compute flow and scaling.
U = E{k}.U1 + E{k}.U2;
V = E{k}.V1 + E{k}.V2;
W = U + V;

G = figure;
axis off;
daspect([1, 1, 1]);
hold on;
colormap(cmap);
set(G, 'renderer', 'opengl');
trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), D.f{1}(Vidx), 'EdgeColor', 'none', 'FaceAlpha', fa);
shading interp;
quiver3(P(:, 1), P(:, 2), P(:, 3), W(Fidx, 1), W(Fidx, 2), W(Fidx, 3), scale, 'k', 'LineWidth', lw);
view(2);
axis(lim);
adjustFigure3;

H = figure;
axis off;
daspect([1, 1, 1]);
hold on;
colormap(cmap);
set(H, 'renderer', 'opengl');
trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), D.f{1}(Vidx), 'EdgeColor', 'none', 'FaceAlpha', fa);
shading interp;
quiver3(P(:, 1), P(:, 2), P(:, 3), U(Fidx, 1), U(Fidx, 2), U(Fidx, 3), scale, 'k', 'LineWidth', lw);
view(2);
axis(lim);
adjustFigure3;

I = figure;
axis off;
daspect([1, 1, 1]);
hold on;
colormap(cmap);
set(I, 'renderer', 'opengl');
trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), D.f{1}(Vidx), 'EdgeColor', 'none', 'FaceAlpha', fa);
shading interp;
quiver3(P(:, 1), P(:, 2), P(:, 3), V(Fidx, 1), V(Fidx, 2), V(Fidx, 3), scale, 'k', 'LineWidth', lw);
view(2);
axis(lim);
adjustFigure3;

J = figure;
axis off;
daspect([1, 1, 1]);
hold on;
colormap(cmap);
set(J, 'renderer', 'opengl');
trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), D.f{2}(Vidx), 'EdgeColor', 'none', 'FaceAlpha', fa);
shading interp;
view(2);
axis(lim);
adjustFigure3;

% Save figures.
savefigure(G, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-1-600dpi.png', filename, k)), '-png', '-r600');
savefigure(G, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-1-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(G, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-1-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(G, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-1-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

savefigure(H, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-2-600dpi.png', filename, k)), '-png', '-r600');
savefigure(H, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-2-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(H, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-2-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(H, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-2-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

savefigure(I, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-3-600dpi.png', filename, k)), '-png', '-r600');
savefigure(I, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-3-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(I, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-3-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(I, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-3-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

savefigure(J, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-4-600dpi.png', filename, k)), '-png', '-r600');
savefigure(J, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-4-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(J, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-4-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(J, fullfile(renderPath, 'detail2', sprintf('%s-%i-detail-4-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

end