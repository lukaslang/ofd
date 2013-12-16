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

% This is used to render and save vector spherical harmonics for
% demonstration purpose.
clear;
close all;
clc;

mkdir(fullfile('./', 'renderings', 'vspharm'));

% Create triangulation of unit sphere.
ref = 3;
[F, V] = sphTriang(ref);
n = size(F, 1);

% Get triangle incenters.
T = TriRep(F, V);
P = T.incenters;

% Set interval of degrees.
N = 2;

% Run through all degrees and all orders.
for k=N
    % Create spherical harmonics.
    [Y1, Y2] = vspharm(k, F, V);
    % Create scalar spherical harmonics for visualisation.
    Ynj = spharm(k, V);
    H = createFigure3;
    for l=1:2*k+1
        cla;
        f = Ynj(:, l);
        trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
        shading interp;
        view(3);
        % Plot vector field.
        quiver3(P(:, 1), P(:, 2), P(:, 3), Y1(:, l, 1), Y1(:, l, 2), Y1(:, l, 3), 1, 'r');
        quiver3(P(:, 1), P(:, 2), P(:, 3), Y2(:, l, 1), Y2(:, l, 2), Y2(:, l, 3), 1, 'b');
        adjustFigure3;
        set(gca, 'ZTick', -1:0.5:1);
        set(gca, 'CLim', [-1, 1]);
        % Save image.
        savefigure(H, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i.png', k, l)));
    end
end
% Save last figure with colorbar.
colorbar;
adjustFigure3;
set(gca, 'ZTick', -1:0.5:1);
cbar = findobj(H, 'tag', 'Colorbar');
set(cbar, 'YTick', -1:0.25:1);
set(cbar, 'TickLength', [.02 .02], 'YColor', [.3 .3 .3]);
savefigure(H, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar.png', k, l)));