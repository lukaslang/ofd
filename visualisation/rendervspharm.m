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
[Faces, Verts] = sphTriang(ref);
n = size(Faces, 1);

% Get triangle incenters.
T = TriRep(Faces, Verts);
P = T.incenters;

% Set interval of degrees.
N = 2;

% Run through all degrees and all orders.
for k=N
    % Create spherical harmonics.
    [Y1, Y2] = vspharm(k, Faces, Verts);
    % Create scalar spherical harmonics for visualisation.
    Ynj = spharm(k, Verts);
    F = createFigure3;
    for l=1:2*k+1
        cla;
        f = Ynj(:, l);
        trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), f);
        shading interp;
        view(3);
        % Plot vector field.
        quiver3(P(:, 1), P(:, 2), P(:, 3), Y1(:, l, 1), Y1(:, l, 2), Y1(:, l, 3), 1, 'r');
        quiver3(P(:, 1), P(:, 2), P(:, 3), Y2(:, l, 1), Y2(:, l, 2), Y2(:, l, 3), 1, 'b');
        adjustFigure3;
        set(gca, 'ZTick', -1:0.5:1);
        set(gca, 'CLim', [-1, 1]);
        % Save image.
        savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-600dpi.png', k, l)), '-png', '-r600');
        savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-1200dpi.png', k, l)), '-png', '-r1200');
        savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-600dpi.jpg', k, l)), '-jpg', '-r600', '-q100');
        savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-1200dpi.jpg', k, l)), '-jpg', '-r1200', '-q100');
    end
end
% Save last figure with colorbar.
colorbar;
adjustFigure3;
set(gca, 'ZTick', -1:0.5:1);
cbar = findobj(F, 'tag', 'Colorbar');
set(cbar, 'YTick', -1:0.25:1);
set(cbar, 'TickLength', [.02 .02], 'YColor', [0 0 0]);
savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-600dpi.png', k, l)), '-png', '-r600');
savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-1200dpi.png', k, l)), '-png', '-r1200');
savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-600dpi.jpg', k, l)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-1200dpi.jpg', k, l)), '-jpg', '-r1200', '-q100');