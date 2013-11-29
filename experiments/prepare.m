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
filename = 'frames-114-116-filtered.mat';
%filename = 'frames-114-116-unfiltered.mat';
%filename = 'frames-56-58-filtered.mat';
%filename = 'frames-56-57-filtered.mat';
% Set working directory.
path = fullfile('./', 'data', name);

% Import data.
disp('Loading image data.');
load(fullfile(path, filename));

% Import cell centres.
disp('Loading cell centres.');
load(fullfile(path, 'thresholdedcenters.mat'));

% Load colormap for proper visualisation.
load(fullfile(path, 'cmapblue.mat'));

% Define cell centres.
frame = 114;
%frame = 58;

% Scaling in z-direction.
zscale = 4.2832;

% Set degrees of bases.
M = 1:5;
N = 6:10;

% Finite difference time parameter.
h = 1;

% Numerical integration tolerance.
tol = 1e-6;

% Spherical band parameters for data.
band = 40;
%band = 60;

% Prepare cell centres.
X = F{frame}.X;
Y = F{frame}.Y;
Z = -zscale * F{frame}.Z;
shift = -min(Z);

% Fit sphere.
[c, r] = sphereFit([X(:), Y(:), Z(:) + shift]);

% Create triangulation of northern hemisphere of the fitted sphere.
ref = 7;
[Faces, Verts] = halfsphTriang(ref);

f = cell(2);
for k=1:2
    % Prepare data.
    u = flipdim(U{k}.u, 3);

    % Project data.
    [um, un, uo] = size(u);
    [X, Y, Z] = ndgrid(1:um, 1:un, 1:uo);

    % Compute radial maximum intensity projection.
    rs = linspace(r-band, r+band, 2*band);
    VB = kron(rs', Verts);
    fb = dataFromCube(c(1)+VB(:, 1), c(2)+VB(:, 2), c(3)+VB(:, 3), X, Y, zscale * Z, u);
    f{k} = max(reshape(fb, size(Verts, 1), length(rs)), [], 2);
end

% Free memory.
clear VB;
clear U;
clear u;
clear X;
clear Y;
clear Z;
clear fb;

% Create output directory.
path = fullfile(path, 'generated');
mkdir(path);

[~, file, ~] = fileparts(filename);
if(all(M == N))
    disp('Compute linear system for same basis, use faster implementation.');
    [dim, U, d, b] = linearsystem(Faces, Verts, N, f{1}, f{2}, h, tol);
    % Define output file.
    genFile = fullfile(path, sprintf('gen-%s-%i-%i-%i.mat', file, N(1), N(end), ref));
    datFile = fullfile(path, sprintf('dat-%s-%i-%i-%i.mat', file, N(1), N(end), ref));
    % Write output.
    disp('Saving generated data.');
    save(datFile, 'Faces', 'Verts', 'f', 'N', 'ref', 'c', 'r', 'h', 'name', 'tol', '-v7.3');
    save(genFile, 'dim', 'U', 'd', 'b', '-v7.3');
else
    disp('Compute linear system for different basis.');
    [dim1, dim2, U, V, W, d1, d2, b] = linearsystemdb(Faces, Verts, M, N, f{1}, f{2}, h, tol);
    % Define output file.
    genFile = fullfile(path, sprintf('gen-%s-%i-%i-%i-%i-%i.mat', file, M(1), M(end), N(1), N(end), ref));
    datFile = fullfile(path, sprintf('dat-%s-%i-%i-%i-%i-%i.mat', file, M(1), M(end), N(1), N(end), ref));
    % Write output.
    disp('Saving generated data.');
    save(datFile, 'Faces', 'Verts', 'f', 'M', 'N', 'ref', 'c', 'r', 'h', 'name', 'tol', '-v7.3');
    save(genFile, 'dim1', 'dim2', 'U', 'V', 'W', 'd1', 'd2', 'b', '-v7.3');
end