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
function U = ofc(N, F, V, f1, f2, h, alpha, s)
%OFC Computes the optical flow on the sphere using the continuity equation.
%
%   U = OFC(N, F, V, f1, f2, h, alpha) takes a triangulation F, V and images
%   f1, f2 on the vertices of the triangulation and returns the optical 
%   flow U assuming the continuity equation. Scalar h is a spacing 
%   parameter and alpha is the regularisation parameter.
%
%   U = OFC(N, F, V, f1, f2, h, alpha, s) takes an additional real scalar s 
%   as the parameter of the Sobolev space H^{s}(S, TS).
%
%   U is defined on the faces F and is of size [size(F, 1), 3].

m = size(V, 1);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(isvector(f1));
assert(isvector(f2));
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);
assert(alpha >= 0);
assert(h > 0);
assert(N > 0);
if(nargin == 8)
    assert(isscalar(s));
else
    s = 1;
end

% Compute linear system based on the continuity equation.
[dim, U, d, b] = linearsystemc(F, V, 1:N, f1, f2, h);

% Solve linear system.
[u, ~] = ofsolve(dim, U, b, d, alpha, s, 30);

% Recover function.
[U1, U2] = vspharmsynth(1:N, F, V, u);

% Add up vector fields of both types.
U = U1 + U2;

end