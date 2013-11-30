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
function plotcoefficients(E, decomp)
%PLOTCOEFFICIENTS Creates barplots of the spherical harmonics coefficients
%of experiments.
%
%   PLOTCOEFFICIENTS(E, decomp) takes a cell array E and creates barplots.
%   If decomp == true then two columns are created, else only one.

for k=1:length(E)
    if(mod(k, 10) == 1)
        figure;
    end
    fprintf('Plotting coefficients %d/%d\n', k, length(E));    
    if(decomp)
        subplot(10, 2, mod(2*(k-1), 20)+1);
        hold on;
        title(sprintf('u, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        bar(E{k}.u);
        subplot(10, 2, mod(2*(k-1), 20)+2);
        hold on;
        title(sprintf('v, s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta));
        bar(E{k}.v);
    else
        subplot(10, 1, mod(k-1, 10)+1);
        hold on;
        title(sprintf('u, s=%g, alpha=%g', E{k}.s, E{k}.alpha));
        bar(E{k}.u);
    end
end

end