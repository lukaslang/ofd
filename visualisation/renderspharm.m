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

% This is used to render and save scalar spherical harmonics for
% demonstration purpose.
clear;
close all;
clc;

mkdir(fullfile('./', 'renderings', 'spharm'));

% Create triangulation of unit sphere.
ref = 5;
[F, V] = sphTriang(ref);

% Set interval of degrees.
N = 0:3;

% Run through all degrees and all orders.
for k=N
    % Create spherical harmonics.
    Ynj = spharm(k, V);
    H = createFigure3;
    for l=1:2*k + 1
        cla;
        trisurf(F, V(:, 1), V(:, 2), V(:, 3), Ynj(:, l), 'EdgeColor', 'none');
        shading interp;
        view(3);
        adjustFigure3;
        set(gca, 'ZTick', -1:0.5:1);
        % Save image.
        savefigure(H, fullfile('./', 'renderings', 'spharm', sprintf('spharm-deg-%i-ord-%i.png', k, l)));
    end
end