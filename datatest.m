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
load(fullfile(path, 'frames-114-116-filtered.mat'));

% Import cell centres.
disp('Loading cell centres.');
load(fullfile(path, 'thresholdedcenters.mat'));

% Load colormap for proper visualisation.
load(fullfile(path, 'cmapblue.mat'));

frame = 114;
k = 1;

% Prepare cell centres.
X = F{frame}.X;
Y = F{frame}.Y;
Z = -4.2832 * F{frame}.Z;
shift = -min(Z);

% Fit sphere.
[c, r] = sphereFit([X(:), Y(:), Z(:) + shift]);

% Prepare data.
u = flipdim(U{k}.u, 3);

% Create parametrisation.
m = 500;
n = 500;
rs = linspace(r-40, r+40, 80);

% Project data.
[um, un, uo] = size(u);
[X, Y, Z] = ndgrid(1:um, 1:un, 1:uo);
f = maximumIntensity(c, m, n, rs, X, Y, 4.2832 * Z, u);

% Create unit sphere for visualisation.
[x, y, z] = sphericalBand(m, n, 1);

% Plot sphere.
figure;
surf(x, y, z, 255 * f, 'Facecolor', 'texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);
figure;
surf(x, y, z, 255 * f, 'Facecolor', 'texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);
colormap(cmap);

% Create parametrisation and polar coordinates.
[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));

% Convert to cartesian coordinates.
r = 1;
x = r .* sqrt(1 - t .^2) .* cos(phi);
y = r .* sqrt(1 - t .^2) .* sin(phi);
z = r .* t;

% Fit spherical harmonics.
N = 10;
[a, Ynj] = spharmFit(N, phi(:), t(:), f(:));

% Recover fitted solution.
g = Ynj*a;
g = reshape(g, m, n);

% Visualise fit.
figure;
surf(x, y, z, g, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

fprintf('Least-squares error is %f.\n', norm(f-g, 2)^2);