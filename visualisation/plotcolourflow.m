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
        % Recover vector field.
        U = E{k}.U1 + E{k}.U2;
        V = E{k}.V1 + E{k}.V2;
        % Compute colour space scaling.
        nmax = max(sqrt(sum((U + V).^2, 2)));

        subplot(1, cols, 3);
        hold on;
        c = double(squeeze(computeColour(U(:, 1)/nmax, U(:, 2)/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('U, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);

        subplot(1, cols, 4);
        hold on;
        c = double(squeeze(computeColour(V(:, 1)/nmax, V(:, 2)/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('V, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);

        subplot(1, cols, 5);
        hold on;
        c = double(squeeze(computeColour((U(:, 1) + V(:, 1))/nmax, (U(:, 2) + V(:, 2))/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('U+V, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);
    else
        % Recover vector field.
        U = E{k}.U1 + E{k}.U2;
        % Compute colour space scaling.
        nmax = max(sqrt(sum(U .^2, 2)));

        subplot(1, cols, 3);
        hold on;
        c = double(squeeze(computeColour(U(:, 1)/nmax, U(:, 2)/nmax))) ./ 255;
        axis([-1, 1, -1, 1, 0, 1]);
        title(sprintf('s=%g, alpha=%g', E{k}.s, E{k}.alpha));
        trisurf(D.Faces, D.Verts(:, 1), D.Verts(:, 2), D.Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
        daspect([1, 1, 1]);
        view(2);
    end
end

end