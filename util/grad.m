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
function gradf = grad(F, V, f, H)
%GRAD Computes the gradient of a function on a triangulation. The function
%is assumed to be piecewise linear and given at the vertices of the
%triangulation.
%
%   gradf = grad(F, V, f) takes a triangulation F, V and a vector f defined
%   on the vertices V and returns the gradient for each triangle in F.
%
%   gradf = grad(F, V, f, H) additionally takes a size(F, 1)-by-3-by-3 
%   matrix H containing the height vectors of the triangles F to be more 
%   efficient in cases where GRAD is called very often.
%
%   Note that size(gradf) = [size(F, 1), 3].

assert(isvector(f));
assert(length(f) == size(V, 1));

if(nargin == 4)
    assert(size(H, 1) == size(F, 1));
end

% Compute height vectors.
if(nargin == 3)
    H = height(F, V);
end
lenH = sum(H(:, 2:3, :).^2, 3);

% Compute values for triangles.
fF = f(F);

% Compute coefficients for height vectors.
c1 = repmat((fF(:, 1) - fF(:, 2)) ./ lenH(:, 1), 1, 3);
c2 = repmat((fF(:, 1) - fF(:, 3)) ./ lenH(:, 2), 1, 3);

% Compute gradient.
gradf = c1 .* squeeze(H(:, 2, :)) + c2 .* squeeze(H(:, 3, :));

end