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
resultsPath = fullfile('./', 'results', name, 'ofd', '2013-11-28-10-10-13');
%resultsPath = fullfile('./', 'results', name, 'ofdb', '2013-11-25-16-18-50');
E = loadexperiments(resultsPath);

% Import data.
disp('Loading precomputed data.');
name = 'cxcr4aMO2_290112';
path = fullfile('./', 'data', name, 'generated');
%filename = 'frames-114-116-unfiltered-1-5-6-10-7';
filename = 'frames-114-116-filtered-1-10-7';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));

% Load colormap for proper visualisation.
load(fullfile('./', 'data', name, 'cmapblue.mat'));

% Plot residual vector.
plotresiduals(E, true);

% Plot coefficients.
plotcoefficients(E, true);

% Plot data and flows.
plotcolourflow(E, D, cmap, true);

% Plot colourwheel.
figure;
cw = colourWheel;
imagesc(cw);