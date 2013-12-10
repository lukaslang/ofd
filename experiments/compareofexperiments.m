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

% Define dataset and get result files.
name = 'cxcr4aMO2_290112';
resultsPath = fullfile('./', 'results', name, 'of');
%load(fullfile(resultsPath, '2013-11-30-10-00-52-frames-114-116-filtered-1-100-7.mat'));
%load(fullfile(resultsPath, '2013-12-03-22-47-45-frames-114-116-filtered-1-100-7.mat'));
%load(fullfile(resultsPath, '2013-12-05-12-08-20-frames-114-116-unfiltered-1-100-7.mat'));
load(fullfile(resultsPath, '2013-12-10-17-02-33-frames-114-116-unfiltered-1-10-7-cont.mat'));

% Import data.
disp('Loading precomputed data.');
name = 'cxcr4aMO2_290112';
path = fullfile('./', 'data', name, 'generated');
%filename = 'frames-114-116-filtered-1-100-7';
%filename = 'frames-114-116-unfiltered-1-100-7';
filename = 'frames-114-116-unfiltered-1-10-7-cont';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));

% Load colormap for proper visualisation.
load(fullfile('./', 'data', name, 'cmapblue.mat'));

% Plot residual vector.
plotresiduals(E, false);

% Plot terms.
plotterms(E, D, false);

% Plot coefficients.
plotcoefficients(E, false);

% Plot data and flows.
plotcolourflow(E, D, cmap, false);

% Plot Helmholtz decomposition.
plothelmholtz(E, D, cmap, false);

% Plot colourwheel.
figure;
cw = colourWheel;
imagesc(cw);