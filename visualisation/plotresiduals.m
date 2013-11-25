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
function plotresiduals(E, decomp)
%PLOTRESIDUALS Plots the residuals of the iterative solves.
%
%   PLOTRESIDUALS(E, decomp) takes a cell vector E and a boolean decomp and
%   plots the residuals.

figure;

% Subplot settings.
cols = min(6, length(E));
rows = ceil(length(E)/cols);

for k=1:length(E)
    fprintf('Plotting residual %d/%d\n', k, length(E));
    subplot(rows, cols, k);
    hold on;
    if(decomp)
        title(sprintf('s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
    else
        title(sprintf('s=%g, alpha=%g', E{k}.s, E{k}.alpha));
    end
    if(strcmp(E{k}.L.solver, 'cgs'))
        plot(0:length(E{k}.L.resvec)-1, E{k}.L.resvec/E{k}.L.rhs, 'b-');
        hold on;
        plot(E{k}.L.iter, E{k}.L.relres, 'rx');
        text(E{k}.L.iter, E{k}.L.relres, sprintf('%0.5f', E{k}.L.relres), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    elseif(strcmp(E{k}.L.solver, 'gmres'))
        plot(0:length(E{k}.L.resvec)-1, E{k}.L.resvec/E{k}.L.rhs, 'b-');
        hold on;
        if(isempty(E{k}.L.restart))
            pos = E{k}.L.iter(2);
        else
            pos = (E{k}.L.iter(1)-1)*E{k}.L.restart+E{k}.L.iter(2);
        end
        plot(pos, E{k}.L.relres, 'rx');
        text(pos, E{k}.L.relres,  sprintf('%0.5f', E{k}.L.relres), 'horizontal', 'right', 'vertical', 'bottom');
    end
end

end