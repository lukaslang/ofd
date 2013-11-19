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
load(fullfile(genPath, 'gen-30-7.mat'));
load(fullfile(genPath, 'dat-30-7.mat'));

% Create folder for results.
resultsPath = fullfile('./', 'results', name, 'of', datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
mkdir(resultsPath);

% Set range for parameters.
rng = [0.01, 0.1, 1, 10, 100, 1000];

% Run experiments.
run = 1;
runs = length(rng);
for alpha=rng
    fprintf('Computing flow %d/%d: %0.2f-cgs\n', run, runs, alpha);
    ticId = tic;
    [U, u, L] = ofsolve(dim, At, b, Y, d, alpha);
    elapsedTime = toc(ticId);
    fprintf('Elapsed time %d seconds.\n', elapsedTime);

    % Create filename.
    wsFilename = sprintf('of-%s-%0.2f-cgs.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), alpha);
    % Save workspace.
    save(fullfile(resultsPath, wsFilename), 'U', 'u', 'L', 'alpha', '-v7.3');
    run = run + 1;
end