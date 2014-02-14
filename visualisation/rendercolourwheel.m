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

% Create render path.
renderPath = fullfile('./', 'renderings', 'colourwheel');
mkdir(renderPath);

% Create colour wheel.
cw = colourWheel;
[m, n, ~] = size(cw);

% Create and normalise coordiantes.
[X, Y, Z] = meshgrid(1:m, 1:n, 0);
X = 2*X / m - 1;
Y = 2*Y / m - 1;

% Get indices.
idx = sqrt(X.^2 + Y.^2) >= 0.96;
[r, c] = ind2sub([m, n], idx);
cwr = cw(:, :, 1);
cwg = cw(:, :, 2);
cwb = cw(:, :, 3);
cwr(idx) = 0;
cwg(idx) = 0;
cwb(idx) = 0;
cw = cat(3, cwr, cwg, cwb);

% Render 3D colour wheel.
F = createFigure3;
set(F, 'Renderer', 'opengl');
surf(X, Y, Z, cw, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'AlphaData', double(~idx), 'AlphaDataMapping', 'none', 'FaceAlpha', 'texturemap');
view(3);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
set(gca, 'ZLim', [0, 1e-3]);
adjustFigure3;
axis off;
file = fullfile(renderPath, 'colourwheel3-600dpi.png');
export_fig(file, '-png', '-r600', '-transparent', F);
file = fullfile(renderPath, 'colourwheel3-1200dpi.png');
export_fig(file, '-png', '-r1200', '-transparent', F);
file = fullfile(renderPath, 'colourwheel3-600dpi.jpg');
export_fig(file, '-jpg', '-r600', '-q100', '-transparent', F);
file = fullfile(renderPath, 'colourwheel3-1200dpi.jpg');
export_fig(file, '-jpg', '-r1200', '-q100', '-transparent', F);

% Save 2D colour wheel.
view(2);
file = fullfile(renderPath, 'colourwheel2-600dpi.png');
export_fig(file, '-png', '-r600', '-transparent', F);
file = fullfile(renderPath, 'colourwheel2-1200dpi.png');
export_fig(file, '-png', '-r1200', '-transparent', F);
file = fullfile(renderPath, 'colourwheel2-600dpi.jpg');
export_fig(file, '-jpg', '-r600', '-q100', '-transparent', F);
file = fullfile(renderPath, 'colourwheel2-1200dpi.jpg');
export_fig(file, '-jpg', '-r1200', '-q100', '-transparent', F);