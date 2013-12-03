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
resultsPath = fullfile('./', 'results', name, 'ofd');
mkdir(resultsPath);

% Specify memory to use.
mem = 45e9;

% Set range for Sobolev parameter s1.
rng1 = [1, 2];
% Set range for Sobolev parameter s2.
rng2 = [-2, -1];
% Set range for alpha.
rng3 = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
% Set range for beta.
rng4 = [1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7];

% Run experiments.
run = 1;
runs = length(rng1)*length(rng2)*length(rng3)*length(rng4);
E = cell(runs, 1);
for s1=rng1
    for s2=rng2
        for alpha=rng3
            for beta=rng4
                fprintf('Computing flow %d/%d: %g-%g-%g-%g-cgs\n', run, runs, s1, s2, alpha, beta);
                ticId = tic;
                [u, v, L] = ofdsolve(G.dim, G.U, G.b, G.d, alpha, beta, s1, s2);
                elapsedTime = toc(ticId);
                fprintf('Elapsed time %d seconds.\n', elapsedTime);

                % Store experiment.
                E{run}.u = u;
                E{run}.v = v;
                E{run}.L = L;
                E{run}.alpha = alpha;
                E{run}.beta = beta;
                E{run}.s1 = s1;
                E{run}.s2 = s2;

                run = run + 1;
            end
        end
    end
end

disp('Recovering vector fields.');
ticId = tic;
e = cell2mat(E);
% Compute vector spherical harmonics synthesis.
[U1, U2] = vspharmsynth(D.N, D.Faces, D.Verts, [e.u, e.v], mem);
elapsedTime = toc(ticId);
fprintf('Elapsed time %d seconds.\n', elapsedTime);

% Save to experiments.
for k=1:runs
    E{k}.U1 = U1(:, :, k);
    E{k}.U2 = U2(:, :, k);
    E{k}.V1 = U1(:, :, k+runs);
    E{k}.V2 = U2(:, :, k+runs);
end

% Create filename.
wsFilename = sprintf('%s-%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), filename);
% Save workspace.
save(fullfile(resultsPath, wsFilename), 'E', '-v7.3');