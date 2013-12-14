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

% This is used to render and save the colour wheel for demonstration.
clear;
close all;
clc;

mkdir(fullfile('./', 'renderings', 'colourwheel'));

% Create colour wheel.
cw = colourWheel;
[m, n, ~] = size(cw);

% Create and normalise coordiantes.
[X, Y] = meshgrid(1:m, 1:n);
X = 2*X / m - 1;
Y = 2*Y / m - 1;

% Render 3D colour wheel.
F = createFigure3;
surf(X, Y, zeros(m, n), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
view(3);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
set(gca, 'ZLim', [0, 1]);
adjustFigure3;
savefigure(F, fullfile('./', 'renderings', 'colourwheel', 'colourwheel3.png'));

% Save 2D colour wheel.
view(2);
savefigure(F, fullfile('./', 'renderings', 'colourwheel', 'colourwheel2.png'));