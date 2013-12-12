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
%filename = 'frames-114-116-filtered-1-100-7';
filename = 'frames-114-116-unfiltered-1-100-7';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));
G = load(fullfile(path, sprintf('gen-%s.mat', filename)));

% Create folder for results.
resultsPath = fullfile('./', 'results', name, 'ofhd');
mkdir(resultsPath);

% Specify memory to use.
mem = 100e9;

% Set maximum number of iterations.
maxit = 100;

% Set parameters for hierarchical decomposition.
%alpha = 1000*pow2(-(0:15));
alpha = 10*ones(5, 1);
% Get number of steps.
nd = length(alpha);
% Set Sobolev space parameter.
%s = ones(nd, 1);
s = [1, 0.5, 0, -0.5, -1];

% Initialise experiments.
E = cell(nd, 1);

% Initialise coefficients matrix.
ud = zeros(G.dim, nd);

% Decompose into components.
b = G.b;
for k=1:nd
    fprintf('Computing flow %d/%d: s=%g, alpha=%g.\n', k, nd, s(k), alpha(k));
    ticId = tic;
    % Solve linear system.
    [u, L] = ofsolve(G.dim, G.U, b, G.d, alpha(k), s(k), maxit);
    elapsedTime = toc(ticId);
    fprintf('Elapsed time %d seconds.\n', elapsedTime);
    % Store coefficients.
    ud(:, k) = u;    
    % Adjust right hand side.
    b = b - G.U*u;
    
    % Store experiment.
    E{k}.L = L;
    E{k}.alpha = alpha(k);
    E{k}.s = s(k);
end

% Compute cummulative sum of coefficients.
ud = cumsum(ud, 2);

disp('Recovering vector fields.');
ticId = tic;
% Compute vector spherical harmonics synthesis.
[U1, U2] = vspharmsynth(D.N, D.Faces, D.Verts, ud, mem);
elapsedTime = toc(ticId);
fprintf('Elapsed time %d seconds.\n', elapsedTime);

% Save to experiments.
for k=1:nd
    E{k}.u = ud(:, k);
    E{k}.U1 = U1(:, :, k);
    E{k}.U2 = U2(:, :, k);
end

% Create filename.
wsFilename = sprintf('%s-%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), filename);
% Save workspace.
save(fullfile(resultsPath, wsFilename), 'E', '-v7.3');