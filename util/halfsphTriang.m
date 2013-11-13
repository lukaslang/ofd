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
function [F, V] = halfsphTriang(n)
%HALFSPHTRIANG Computes a triangulation of northern half of the unit sphere.
%
%   [F, V] = HALFSPHTRIANG returns faces F and vertices V of a 
%   triangulation of the northern half of the unit sphere, i.e. having
%   z-coordinate >= 0.
%
%   [F, V] = HALFSPHTRIANG(n) does n > 0 refinements of the mesh.
%
%   Note that, for quality reasons of the mesh, if n > 0 the whole unit 
%   sphere is meshed first!
%
%   For details on F, V see sphTriang.

if(nargin == 0)
    n = 0;
end

% Generate unit sphere mesh.
[F, V] = sphTriang(n);

% Get faces on the northern hemisphere.
idx = V(F(:, 1), 3) >= 0 & V(F(:, 2), 3) >= 0 & V(F(:, 3), 3) >= 0;
F = F(idx, :);

% Find vertices on the northern hemisphere.
ids = find(V(:, 3) >= 0);
V = V(ids, :);

% Compute new vertex indices.
nids(ids) = 1:length(ids);
F = nids(F);

end