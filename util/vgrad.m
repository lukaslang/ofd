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
function G = vgrad(F, V, Y, Z, H, lenH, FN)
%VGRAD Computes an approximation of the covariant derivative on a
%triangulation.
%
%   G = VGRAD(F, V, Y, Z) takes a triangulation F, V and a vector field Y 
%   defined on the vertices V and returns the covariant derivative of Y 
%   with respect to another vector field Z on each triangle. Note that Z
%   must be tangential and piecewise constant on each face!
%
%   G = VGRAD(F, V, Y, Z, H) additionally takes a matrix H containing the
%   height vectors of the triangles.
%
%   G = VGRAD(F, V, Y, Z, H, lenH) additionally takes a matrix lenH 
%   containing the length height vectors H(:, 2:3, :).
%
%   G = VGRAD(F, V, Y, Z, H, lenH, FN) additionally takes a matrix FN 
%   containing the (outward pointing) face normals.
%
%   The last three are more efficient in cases where VGRAD is called 
%   multiple times on the same triangulation.
%
%   Note that H must be of size size(F, 1)-by-3-by-3, lenH of size
%   size(F, 1)-by-2, and FN of size size(F, 1)-by-3.
%
%   G is of size m-by-3, where m is the number of faces in F.

% Compute height vectors of triangles spanned by incenters.
if(nargin < 5)
    H = height(F, V);
end
% Compute squared length of 2nd and 3rd height vectors.
if(nargin < 6)
    lenH = sum(H(:, 2:3, :).^2, 3);
end
% Compute face normals.
if(nargin < 7)
    T = TriRep(F, V);
    FN = -T.faceNormals;
end

assert(size(Y, 1) == size(V, 1));
assert(size(Y, 2) == 3);
assert(size(Z, 1) == size(F, 1));
assert(size(Z, 2) == 3);
assert(size(H, 1) == size(F, 1));
assert(size(H, 2) == 3);
assert(size(H, 3) == 3);
assert(size(lenH, 1) == size(F, 1));
assert(size(lenH, 2) == 2); 
assert(size(FN, 1) == size(F, 1));
assert(size(FN, 2) == 3);

% Compute componentwise directional derivative D_Z Y^i.
G = zeros(size(F, 1), 3);
parfor k=1:3
    g = grad(F, V, Y(:, k), H, lenH);
    G(:, k) = dot(Z, g, 2); 
end

% Compute normal part.
ip = dot(G, FN, 2);

% Project onto tangent space.
G = G - bsxfun(@times, FN, ip);

end