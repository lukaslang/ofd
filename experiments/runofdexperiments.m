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
filename = 'frames-114-116-filtered-1-10-7';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));
G = load(fullfile(path, sprintf('gen-%s.mat', filename)));

% Create folder for results.
resultsPath = fullfile('./', 'results', name, 'ofd', datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
mkdir(resultsPath);

% Set range for Sobolev parameter s1.
rng1 = [1];
% Set range for Sobolev parameter s2.
rng2 = [-1];
% Set range for alpha.
rng3 = [0.001, 0.01, 0.1, 1];
% Set range for beta.
rng4 = [1, 10, 100, 1000];

% Create vector spherical harmonics.
disp('Generating vector spherical harmonics.');
ticId = tic;
[Y, ~] = vspharmn(D.N, D.Faces, D.Verts);
elapsedTime = toc(ticId);
fprintf('Elapsed time %d seconds.\n', elapsedTime);

% Run experiments.
run = 1;
runs = length(rng1)*length(rng2)*length(rng3)*length(rng4);
for s1=rng1
    for s2=rng2
        for alpha=rng3
            for beta=rng4
                fprintf('Computing flow %d/%d: %g-%g-%g-%g-cgs\n', run, runs, s1, s2, alpha, beta);
                ticId = tic;
                [u, v, L] = ofdsolve(G.dim, G.U, G.b, G.d, alpha, beta, s1, s2);
                elapsedTime = toc(ticId);
                fprintf('Elapsed time %d seconds.\n', elapsedTime);

                disp('Recovering vector field.');
                ticId = tic;
                U = synth(Y, u);
                V = synth(Y, v);
                elapsedTime = toc(ticId);
                fprintf('Elapsed time %d seconds.\n', elapsedTime);
                
                % Create filename.
                wsFilename = sprintf('%s-%s-%g-%g-%g-%g-%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), filename, s1, s2, alpha, beta, L.solver);
                % Save workspace.
                save(fullfile(resultsPath, wsFilename), 'U', 'V', 'u', 'v', 'L', 'alpha', 'beta', 's1', 's2', '-v7.3');
                run = run + 1;
            end
        end
    end
end