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
function [nF, nV] = refineMesh(F, V)
%REFINEMESH Refines the given mesh.
%
%   [F, V] = refineMesh(F, V) takes triangular faces F and vertices V and 
%   refines the mesh by dividing each triangle into four triangles. F is 
%   assumed to contain clockwise indexed triangles and the resulting 
%   triangles are also oriented clockwise.
%
%   Note that the newly inserted vertices still lie on the edges of the 
%   old triangles and must be adjusted to the geometry!

assert(~isempty(F));
assert(~isempty(V));
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);

numF = size(F, 1);
numV = size(V, 1);

% Create adjacency matrix.
TR = TriRep(F, V);
E = TR.edges;
clear TR;

% Assumes vertices in ascending order.
M = sparse(E(:, 1), E(:, 2), numV+1:numV+size(E, 1), numV, numV);

% Initialise new faces.
nF = zeros(4*numF, 3);
% Initialise coordinates of new vertices.
nV = [V; zeros(size(E, 1), 3)];

% Iterate over all faces.
for k=1:numF
    v = V(F(k, :), :);
    % Compute vertex indices of midpoints.   
    idx = [M(min(F(k, 1), F(k, 2)), max(F(k, 1), F(k, 2))); 
        M(min(F(k, 1), F(k, 3)), max(F(k, 1), F(k, 3))); 
        M(min(F(k, 2), F(k, 3)), max(F(k, 2), F(k, 3)))];
    % Compute coordinates of midpoints.
    nV(idx, :) = [v(1, :) + v(2, :); v(1, :) + v(3, :); v(2, :) + v(3, :)]/2;
       
    % Create new faces.
    nF(4*k-3, :) = [F(k, 1), idx(1), idx(2)];
    nF(4*k-2, :) = [idx(2), idx(1), idx(3)];
    nF(4*k-1, :) = [F(k, 3), idx(2), idx(3)];
    nF(4*k, :) = [idx(3), idx(1), F(k, 2)];
end
end