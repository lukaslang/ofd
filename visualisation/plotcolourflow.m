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
function plotcolourflow(E, D, cmap, decomp)
%PLOTCOLOURFLOW Plots the flow of experiments.
%
%   PLOTCOLOURFLOW(E, D, cmap, decomp) takes a cell array E, data D, a colourmap
%   cmap. If decomp == true then five columns are created, else only three.

if(decomp)
    cols = 5;
else
    cols = 3;
end

for k=1:length(E)
    figure;
    colormap(cmap);
    for l=1:2
        subplot(1, cols, l);
        hold on;
        axis([-1, 1, -1, 1, 0, 1]);
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), D.f{l}, 'EdgeColor', 'none');
        shading interp;
        daspect([1, 1, 1]);
        view(2);
    end
    fprintf('Plotting flow %d/%d\n', k, length(E));
    if(decomp)
        % Compute colour space scaling.
        nmax = max(sqrt(sum((E{k}.U + E{k}.V).^2, 2)));

        subplot(1, cols, 3);
        hold on;
        c = double(squeeze(computeColour(E{k}.U(:, 1)/nmax, E{k}.U(:, 2)/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('U, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);

        subplot(1, cols, 4);
        hold on;
        c = double(squeeze(computeColour(E{k}.V(:, 1)/nmax, E{k}.V(:, 2)/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('V, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);

        subplot(1, cols, 5);
        hold on;
        c = double(squeeze(computeColour((E{k}.U(:, 1) + E{k}.V(:, 1))/nmax, (E{k}.U(:, 2) + E{k}.V(:, 2))/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('U+V, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);
    else
        % Compute colour space scaling.
        nmax = max(sqrt(sum(E{k}.U .^2, 2)));

        subplot(1, cols, 3);
        hold on;
        c = double(squeeze(computeColour(E{k}.U(:, 1)/nmax, E{k}.U(:, 2)/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('s=%g, alpha=%g', E{k}.s, E{k}.alpha));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);
    end
end

end