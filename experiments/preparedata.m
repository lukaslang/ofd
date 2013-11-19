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

% Set parameters.
N = 30;
h = 1;
tol = 1e-6;

% Prepare cell centres.
X = F{frame}.X;
Y = F{frame}.Y;
Z = -4.2832 * F{frame}.Z;
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
    rs = linspace(r-40, r+40, 80);
    VB = kron(rs', Verts);
    fb = dataFromCube(c(1)+VB(:, 1), c(2)+VB(:, 2), c(3)+VB(:, 3), X, Y, 4.2832 * Z, u);
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

% Compute data function.
disp('Compute data functions.');
[dim, At, d, Y, b] = computeDataFunctions(Faces, Verts, N, f{1}, f{2}, h, tol);

% Create directory.
path = fullfile(path, 'generated');
mkdir(path);

% Define output files.
dataOut = sprintf('dat-%i-%i.mat', N, ref);
genOut = sprintf('gen-%i-%i.mat', N, ref);

% Write output.
disp('Saving generated data.');
save(fullfile(path, dataOut), 'Faces', 'Verts', 'f', 'N', 'ref', 'c', 'r', 'h', 'name', 'tol', '-v7.3');
save(fullfile(path, genOut), 'dim', 'At', 'd', 'Y', 'b', '-v7.3');