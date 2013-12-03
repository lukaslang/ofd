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
function [U, V] = ofdb(M, N, F, V, f1, f2, h, alpha, beta, s1, s2)
%OFDB Returns an optical flow decomposition on a triangulation using
%different bases for U and V.
%
%   [U, V] = OFDB(M, N, F, V, f1, f2, h, alpha, beta) takes a triangulation 
%   F, V and images f1, f2 on the vertices of the triangulation and returns
%   an optical flow decomposition U, V. Scalar h is a spacing parameter and
%   alpha, beta are regularisation parameters.
%
%   [U, V] = OFDB(M, N, F, V, f1, f2, h, alpha, beta, s1, s2) takes 
%   additional real scalars s1, s2 as the parameters of the Sobolev spaces 
%   H^{s1}(S, TS), H^{s2}(S, TS).
%
%   Degrees M, N must be vectors of consecutive positive integers!
%
%   U and V are vector fields defined on the faces F and are of size 
%   [n, 3], where n = size(F, 1) is the number of faces.
%
%   In case M equals N, see ofd for a more efficient implementation!

m = size(V, 1);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(isvector(f1));
assert(isvector(f2));
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);
assert(alpha >= 0);
assert(beta >= 0);
assert(h > 0);
if(nargin == 11)
    assert(isscalar(s1));
    assert(isscalar(s2));
elseif(nargin == 9)
    s1 = 1;
    s2 = -1;
end

% Compute linear system.
[dim1, dim2, Um, Vm, Wm, d1, d2, b] = linearsystemdb(F, V, M, N, f1, f2, h, 1e-6);

% Solve linear system.
[u, v, ~] = ofdbsolve(dim1, dim2, Um, Vm, Wm, d1, d2, b, alpha, beta, s1, s2, 30);

% Recover functions.
[U1, U2] = vspharmsynth(M, F, V, u);
[V1, V2] = vspharmsynth(N, F, V, v);

% Add up vector fields of both types.
U = U1 + U2;
V = V1 + V2;

end