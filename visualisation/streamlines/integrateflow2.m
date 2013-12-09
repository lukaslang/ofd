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
function verts = integrateflow2(X, V, S, h, maxit)
%INTEGRATEFLOW2 Integrates a static velocity field in the unit circle.
%
%   verts = INTEGRATEFLOW2(X, V, S, h, maxit) takes a velocity field V at 
%   points X and initial values S and returns a cell array verts, where 
%   verts{k} contains the coordinates according to S(k, :). Scalar h is the
%   step size and maxit is the maximum number of steps.

assert(size(X, 2) == 2);
assert(size(V, 2) == 2);
assert(size(S, 2) == 2);
assert(size(X, 1) == size(V, 1));
assert(maxit > 0);
assert(h > 0);

% Initialise vertices.
verts = cell(size(S, 1), 1);

% Interpolate vector field.
F1 = TriScatteredInterp(X, V(:, 1));
F2 = TriScatteredInterp(X, V(:, 2));

for k=1:length(verts)
    posx = S(k, 1);
    posy = S(k, 2);
    verts{k} = [posx, posy];
    for l=1:maxit
        x = posx + h * F1(posx, posy);
        y = posy + h * F2(posx, posy);
        % Stop if trajectory leaves unit circle.
        if(x^2 + y^2 <= 1)        
            verts{k} = [verts{k}; x, y];
            posx = x;
            posy = y;
        else
            break;
        end
    end
end

end