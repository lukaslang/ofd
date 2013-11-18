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
function [U, V] = ofd(N, F, V, f1, f2, h, alpha, beta)
%OFD Returns an optical flow decomposition on a triangulation.
%
%   [U, V] = OFD(N, F, V, f1, f2, h, alpha, beta) takes a triangulation 
%   F, V and images f1, f2 on the vertices of the triangulation and returns
%   an optical flow decomposition U, V. Scalar h is a spacing parameter and
%   alpha, beta are regularisation parameters.
%
%   U and V are defined on the faces F and are of size [n, 3], where 
%   n = size(F, 1) is the number of faces.

m = size(V, 1);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(isvector(f1));
assert(isvector(f2));
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);
assert(alpha > 0);
assert(beta > 0);
assert(h > 0);

% Compute functions needed for solving the linear system.
[dim, At, d, Y, b] = computeDataFunctions(F, V, N, f1, f2, h, 1e-6);

% Solve linear system.
[U, V, ~, ~] = ofdsolve(dim, At, b, Y, d, alpha, beta);

end