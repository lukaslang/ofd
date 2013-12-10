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
function plotterms(E, D, decomp)
%PLOTTERMS Plots the optical flow and regularisation term residuals.
%
%   PLOTTERMS(E, D, decomp) takes a cell vector E, data D, and a boolean 
%   decomp and plots the terms of the functional.

% Compute terms.
r = zeros(length(E), 1);
m = zeros(length(E), 1);
n = zeros(length(E), 1);
o = zeros(length(E), 1);
labels = cell(length(E), 1);
for k=1:length(E)
    fprintf('Computing terms %i/%i.\n', k, length(E));
    % Compute optical flow residual.
    [res, min] = residual(E{k}.U1 + E{k}.U2, D.Faces, D.Verts, D.f{1}, D.f{2}, D.tol);
    r(k) = res;
    m(k) = min;
    % Compute regularisation terms.
    if(decomp)
        n(k) = E{k}.alpha * hnorm(E{k}.s1, vspharmeigs(D.N), E{k}.u);
        o(k) = E{k}.beta * hnorm(E{k}.s2, vspharmeigs(D.N), E{k}.v);
        labels{k} = sprintf('s1=%g, s2=%g, alpha=%g, beta=%g', E{k}.s1, E{k}.s2, E{k}.alpha, E{k}.beta);
    else
        n(k) = E{k}.alpha * hnorm(E{k}.s, vspharmeigs(D.N), E{k}.u);
        labels{k} = sprintf('s=%g, alpha=%g', E{k}.s, E{k}.alpha);
    end
end

% Determine label positions.
lbx = 1:length(E);
figure;
hold on;
plot(lbx, r, '-*r');
plot(lbx, m, '--*b');
plot(lbx, n, '-.*g');
if(decomp)
    plot(lbx, o, '-.*m');
    lby = max([r, m ,n, o], [], 2);
else
    lby = max([r, m ,n], [], 2);
end
text(lbx, lby, labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Rotation', 90, 'FontSize', 8);

end