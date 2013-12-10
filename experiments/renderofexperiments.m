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
resultsPath = fullfile('./', 'results', name, 'of');
%resultsname = '2013-12-03-19-32-46-frames-114-116-filtered-1-10-7';
%resultsname = '2013-12-05-12-08-20-frames-114-116-unfiltered-1-100-7';
resultsname = '2013-12-10-17-02-33-frames-114-116-unfiltered-1-10-7-cont';
load(fullfile(resultsPath, sprintf('%s.mat', resultsname)));

% Import data.
disp('Loading precomputed data.');
name = 'cxcr4aMO2_290112';
path = fullfile('./', 'data', name, 'generated');
%filename = 'frames-114-116-filtered-1-10-7';
filename = 'frames-114-116-unfiltered-1-100-7';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));

% Load colormap for proper visualisation.
load(fullfile('./', 'data', name, 'cmapblue.mat'));

% Define renderings path.
renderPath = fullfile('./', 'renderings', name, 'of', resultsname);
mkdir(renderPath);
mkdir(fullfile(renderPath, 'residual'));
mkdir(fullfile(renderPath, 'coefficients'));
mkdir(fullfile(renderPath, 'data2'));
mkdir(fullfile(renderPath, 'data3'));
mkdir(fullfile(renderPath, 'flow2'));
mkdir(fullfile(renderPath, 'flow3'));
mkdir(fullfile(renderPath, 'helmu2'));
mkdir(fullfile(renderPath, 'helmu3'));
mkdir(fullfile(renderPath, 'helmv2'));
mkdir(fullfile(renderPath, 'helmv3'));

% Restriction allows to search among the results.
%e = cell2mat(E);
%idx = find([e.s] == 1);
%E = E(idx);

for k=1:length(E)

% Plot residual vector.
R = createFigure;
plot(0:length(E{k}.L.resvec)-1, E{k}.L.resvec/E{k}.L.rhs, 'b-');
if(isempty(E{k}.L.restart))
    pos = E{k}.L.iter(2);
else
    pos = (E{k}.L.iter(1)-1)*E{k}.L.restart+E{k}.L.iter(2);
end
plot(pos, E{k}.L.relres, 'rx');
text(pos, E{k}.L.relres,  sprintf('%0.5f', E{k}.L.relres), 'horizontal', 'right', 'vertical', 'bottom');
axis on;
adjustFigure;
savefigure(R, fullfile(renderPath, 'residual', sprintf('%s-%i.png', filename, k)));

% Plot coefficients.
C = createFigure;
bar(E{k}.u);
axis on;
adjustFigure;
savefigure(C, fullfile(renderPath, 'coefficients', sprintf('%s-%i.png', filename, k)));

% Plot data.
for l=1:2
    F = createFigure3(cmap, [-1, 1, -1, 1, 0, 1]);
    trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), D.f{l}, 'EdgeColor', 'none');
    shading interp;
    view(3);
    adjustFigure3;
    savefigure(F, fullfile(renderPath, 'data3', sprintf('%s-%i-%i.png', filename, k, l)));
    
    view(2);
    adjustFigure3;
    savefigure(F, fullfile(renderPath, 'data2', sprintf('%s-%i-%i.png', filename, k, l)));
end

% Plot data and flows.
% Recover vector field.
U1 = projecttoplane(E{k}.U1);
U2 = projecttoplane(E{k}.U2);
U = U1 + U2;

% Compute colour space scaling.
nmax = max(sqrt(sum(U.^2, 2)));

V = createFigure3(cmap, [-1, 1, -1, 1, 0, 1]);
c = double(squeeze(computeColour(U(:, 1)/nmax, U(:, 2)/nmax))) ./ 255;
trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
view(3);
adjustFigure3;
savefigure(V, fullfile(renderPath, 'flow3', sprintf('%s-%i.png', filename, k)));
view(2);
adjustFigure3;
savefigure(V, fullfile(renderPath, 'flow2', sprintf('%s-%i.png', filename, k)));

% Plot Helmholtz decomposition.
H = createFigure3(cmap, [-1, 1, -1, 1, 0, 1]);
c = double(squeeze(computeColour(U1(:, 1)/nmax, U1(:, 2)/nmax))) ./ 255;
trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
view(3);
adjustFigure3;
savefigure(H, fullfile(renderPath, 'helmu3', sprintf('%s-%i.png', filename, k)));
view(2);
adjustFigure3;
savefigure(H, fullfile(renderPath, 'helmu2', sprintf('%s-%i.png', filename, k)));

I = createFigure3(cmap, [-1, 1, -1, 1, 0, 1]);
c = double(squeeze(computeColour(U2(:, 1)/nmax, U2(:, 2)/nmax))) ./ 255;
trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
view(3);
adjustFigure3;
savefigure(I, fullfile(renderPath, 'helmv3', sprintf('%s-%i.png', filename, k)));
view(2);
adjustFigure3;
savefigure(I, fullfile(renderPath, 'helmv2', sprintf('%s-%i.png', filename, k)));

close all;
end