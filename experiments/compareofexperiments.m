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
resultsPath = fullfile('./', 'results', name, 'of', '2013-11-19-10-16-14');
files = getFiles(resultsPath);

% Import data.
disp('Loading precomputed data.');
name = 'cxcr4aMO2_290112';
genPath = fullfile('./', 'data', name, 'generated');
load(fullfile(genPath, 'dat-30-7.mat'));

% Load colormap for proper visualisation.
load(fullfile('./', 'data', name, 'cmapblue.mat'));

% Load experiments.
E = cell(length(files), 1);
for k=1:length(files)
    filename = files(k).name;
    fprintf('Loading run %d/%d: %s\n', k, length(files), filename);
    E{k} = load(fullfile(resultsPath, filename));
end

% Subplot settings.
cols = min(6, length(files));
rows = ceil(length(files)/cols);

% Plot residual vector.
figure;
for k=1:length(files)
    fprintf('Plotting residual %d/%d\n', k, length(files));
    subplot(rows, cols, k);
    hold on;
    title(sprintf('alpha=%g', E{k}.alpha));
    if(strcmp(E{k}.L.solver, 'cgs'))
        plot(0:length(E{k}.L.resvec)-1, E{k}.L.resvec/E{k}.L.rhs, 'b-');
        hold on;
        plot(E{k}.L.iter, E{k}.L.relres, 'rx');
        text(E{k}.L.iter, E{k}.L.relres, sprintf('%0.5f', E{k}.L.relres), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    elseif(strcmp(E{k}.L.solver, 'gmres'))
        plot(0:length(E{k}.L.resvec)-1, E{k}.L.resvec/E{k}.L.rhs, 'b-');
        hold on;
        if(isempty(E{k}.L.restart))
            pos = E{k}.L.iter(2);
        else
            pos = (E{k}.L.iter(1)-1)*E{k}.L.restart+E{k}.L.iter(2);
        end
        plot(pos, E{k}.L.relres, 'rx');
        text(pos, E{k}.L.relres,  sprintf('%0.5f', E{k}.L.relres), 'horizontal', 'right', 'vertical', 'bottom');
    end
end

% Plot coefficients.
figure;
for k=1:length(files)
    fprintf('Plotting coefficients %d/%d\n', k, length(files));
    subplot(length(files), 1, k);
    hold on;
    title(sprintf('u, alpha=%g', E{k}.alpha));
    bar(E{k}.u);
end

% Plot data and flow.
for k=1:length(files)
    figure;
    colormap(cmap);
    for l=1:2
        subplot(1, 3, l);
        hold on;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('alpha=%g', E{k}.alpha));
        trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), f{l}, 'EdgeColor', 'none');
        shading interp;
        daspect([1, 1, 1]);
        view(2);
    end
    fprintf('Plotting flow %d/%d\n', k, length(files));
    % Compute colour space scaling.
    nmax = max(sqrt(sum(E{k}.U.^2, 2)));
    % Compute colour of projection.
    c = double(squeeze(computeColour(E{k}.U(:, 1)/nmax, E{k}.U(:, 2)/nmax))) ./ 255;
    subplot(1, 3, 3);
    hold on;
    axis([-1, 1, -1, 1, 0, 1]);
    title(sprintf('alpha=%g', E{k}.alpha));
    trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(2);
end

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(2);