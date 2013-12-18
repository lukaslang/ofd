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
% Optical flow.
%load(fullfile(resultsPath, '2013-12-03-22-47-45-frames-114-116-filtered-1-100-7.mat'));
load(fullfile(resultsPath, '2013-12-05-12-08-20-frames-114-116-unfiltered-1-100-7.mat'));

% Optical flow with continuity equation.
%load(fullfile(resultsPath, '2013-12-10-23-40-35-frames-114-116-unfiltered-1-100-7-cont.mat'));

% Optical flow hierarchical decomposition.
%resultsPath = fullfile('./', 'results', name, 'ofhd');
%load(fullfile(resultsPath, '2013-12-11-21-53-46-frames-114-116-unfiltered-1-100-7.mat'));
%load(fullfile(resultsPath, '2013-12-15-14-55-13-frames-114-116-unfiltered-1-100-7.mat'));

% Import data.
disp('Loading precomputed data.');
path = fullfile('./', 'data', name, 'generated');
%filename = 'frames-114-116-filtered-1-100-7';
filename = 'frames-114-116-unfiltered-1-100-7';
%filename = 'frames-114-116-unfiltered-1-100-7-cont';
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