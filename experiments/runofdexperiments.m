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
genPath = fullfile('./', 'data', name, 'generated');
load(fullfile(genPath, 'gen-50-7.mat'));
load(fullfile(genPath, 'dat-50-7.mat'));

% Create folder for results.
resultsPath = fullfile('./', 'results', name, 'ofd', datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
mkdir(resultsPath);

% Set range for parameters.
rng1 = [0.001, 0.01, 0.1];
rng2 = [1, 10, 100, 1000];

% Run experiments.
run = 1;
runs = length(rng1)*length(rng2);
for alpha=rng1
    for beta=rng2
        fprintf('Computing decomposition %d/%d: %g-%g-cgs\n', run, runs, alpha, beta);
        ticId = tic;
        [U, V, u, v, L] = ofdsolve(dim, At, b, Y, d, alpha, beta);
        elapsedTime = toc(ticId);
        fprintf('Elapsed time %d seconds.\n', elapsedTime);

        % Create filename.
        wsFilename = sprintf('ofd-%s-%g-%g-%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), alpha, beta, L.solver);
        % Save workspace.
        save(fullfile(resultsPath, wsFilename), 'U', 'V', 'u', 'v', 'L', 'alpha', 'beta', '-v7.3');
        run = run + 1;
    end
end