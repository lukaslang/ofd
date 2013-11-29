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
resultsPath = fullfile('./', 'results', name, 'of', datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
mkdir(resultsPath);

% Set range for Sobolev parameter s.
rng1 = [-2, -1, -0.5, 0, 0.5, 1, 2];
% Set range for alpha.
rng2 = [0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000];

% Create vector spherical harmonics.
disp('Generating vector spherical harmonics.');
ticId = tic;
[Y, ~] = vspharmn(D.N, D.Faces, D.Verts);
elapsedTime = toc(ticId);
fprintf('Elapsed time %d seconds.\n', elapsedTime);

% Run experiments.
run = 1;
runs = length(rng1)*length(rng2);
for s=rng1
    for alpha=rng2
        fprintf('Computing flow %d/%d: %g-%g-cgs\n', run, runs, s, alpha);
        ticId = tic;
        [u, L] = ofsolve(G.dim, G.U, G.b, G.d, alpha, s);
        elapsedTime = toc(ticId);
        fprintf('Elapsed time %d seconds.\n', elapsedTime);
        
        disp('Recovering vector field.');
        ticId = tic;
        U = synth(Y, u);
        elapsedTime = toc(ticId);
        fprintf('Elapsed time %d seconds.\n', elapsedTime);
        
        % Create filename.
        wsFilename = sprintf('%s-%s-%g-%g-%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), filename, s, alpha, L.solver);
        % Save workspace.
        save(fullfile(resultsPath, wsFilename), 'U', 'u', 'L', 'alpha', 's', '-v7.3');
        run = run + 1;
    end
end