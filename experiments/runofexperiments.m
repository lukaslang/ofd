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

% Import data.
disp('Loading precomputed data.');
name = 'cxcr4aMO2_290112';
path = fullfile('./', 'data', name, 'generated');
filename = 'frames-114-116-filtered-1-100-7';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));
G = load(fullfile(path, sprintf('gen-%s.mat', filename)));

% Create folder for results.
resultsPath = fullfile('./', 'results', name, 'of');
mkdir(resultsPath);

% Specify memory to use.
mem = 50e9;

% Set maximum number of iterations.
maxit = 100;

% Set range for Sobolev parameter s.
rng1 = [-2, -1, -0.5, 0, 0.5, 1, 2];
% Set range for alpha.
rng2 = [0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000];

% Run experiments.
run = 1;
runs = length(rng1)*length(rng2);
E = cell(runs, 1);
for s=rng1
    for alpha=rng2
        fprintf('Computing flow %d/%d: s=%g, alpha=%g\n', run, runs, s, alpha);
        ticId = tic;
        [u, L] = ofsolve(G.dim, G.U, G.b, G.d, alpha, s, maxit);
        elapsedTime = toc(ticId);
        fprintf('Elapsed time %d seconds.\n', elapsedTime);

        % Store experiment.
        E{run}.u = u;
        E{run}.L = L;
        E{run}.alpha = alpha;
        E{run}.s = s;
        
        run = run + 1;
    end
end

disp('Recovering vector fields.');
ticId = tic;
e = cell2mat(E);
% Compute vector spherical harmonics synthesis.
[U1, U2] = vspharmsynth(D.N, D.Faces, D.Verts, [e.u], mem);
elapsedTime = toc(ticId);
fprintf('Elapsed time %d seconds.\n', elapsedTime);

% Save to experiments.
for k=1:runs
    E{k}.U1 = U1(:, :, k);
    E{k}.U2 = U2(:, :, k);
end

% Create filename.
wsFilename = sprintf('%s-%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), filename);
% Save workspace.
save(fullfile(resultsPath, wsFilename), 'E', '-v7.3');