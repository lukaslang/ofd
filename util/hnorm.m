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
function n = hnorm(s, a, u)
%HNORM Computes the Sobolev norm of the space H^s(a).
%
%   n = HNORM(s, a, u) takes a scalar s and a vector a, which is a sequence,
%   and Fourier coefficients u and returns the norm of a vector field in
%   the space H^s(a).
%
%   Note that a and u must be of equal length. Note that n is a scalar.

assert(isscalar(s));
assert(isvector(a));
assert(isvector(u));
assert(length(a) == length(u));

% Compute induced Sobolev norm using Parseval's identity.
n = sqrt((a'.^s) * (u.^2));

end