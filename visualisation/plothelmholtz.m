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
function plothelmholtz(E, D, cmap, decomp)
%PLOTHELMHOLTZ Plots the Helmholtz decomposition of experiments.
%
%   PLOTHELMHOLTZ(E, D, cmap, decomp) takes a cell array E, data D, a colourmap
%   cmap. If decomp == true then nine columns are created, else only five.

if(decomp)
    cols = 9;
else
    cols = 5;
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
    fprintf('Plotting Helmholtz decomposition %d/%d\n', k, length(E));
    if(decomp)
        % Recover vector field.
        U1 = projecttoplane(E{k}.U1);
        U2 = projecttoplane(E{k}.U2);
        V1 = projecttoplane(E{k}.V1);
        V2 = projecttoplane(E{k}.V2);
        U = U1 + U2;
        V = V1 + V2;
        W = U + V;
        
        % Compute colour space scaling.
        nmax = max(sqrt(sum(W.^2, 2)));

        subplot(1, cols, 3);
        titlestr = sprintf('U1, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
        plot(U1/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 4);
        titlestr = sprintf('U2, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
        plot(U2/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 5);
        titlestr = sprintf('U, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
        plot(U/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 6);
        titlestr = sprintf('V1, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
        plot(V1/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 7);
        titlestr = sprintf('V2, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
        plot(V2/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 8);
        titlestr = sprintf('V, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
        plot(V/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 9);
        titlestr = sprintf('U+V, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
        plot(W/nmax, D.Faces, D.Verts, titlestr);
    else
        % Recover vector field.
        U1 = projecttoplane(E{k}.U1);
        U2 = projecttoplane(E{k}.U2);
        U = U1 + U2;
        
        % Compute colour space scaling.
        nmax = max(sqrt(sum(U .^2, 2)));

        subplot(1, cols, 3);
        titlestr = sprintf('U1, s=%g, alpha=%g', E{k}.s, E{k}.alpha);
        plot(U1/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 4);
        titlestr = sprintf('U2, s=%g, alpha=%g', E{k}.s, E{k}.alpha);
        plot(U2/nmax, D.Faces, D.Verts, titlestr);
        
        subplot(1, cols, 5);
        titlestr = sprintf('U1+U2, s=%g, alpha=%g', E{k}.s, E{k}.alpha);
        plot(U/nmax, D.Faces, D.Verts, titlestr);
    end
end

end

function p = projecttoplane(v)
    p = v;
    p(:, 3) = 0;
    lenv = sqrt(sum(v .^2, 2));
    lenp = sqrt(sum(p .^2, 2));
    p = bsxfun(@times, p, lenv./ lenp);
end

function plot(v, Faces, Verts, titlestr)
    hold on;
    c = double(squeeze(computeColour(v(:, 1), v(:, 2)))) ./ 255;
    axis([-1, 1, -1, 1, 0, 1]);
    title(titlestr);
    trisurf(Faces, Verts(:, 1), Verts(:, 2), Verts(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(2);
end