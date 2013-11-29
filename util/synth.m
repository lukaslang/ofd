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
function U = synth(Y, u)
%SYNTH Takes vector spherical harmonics and coefficients and computes the 
%synthesis.
%
%   U = SYNTH(Y, u) takes vector spherical harmonics Y and coefficients u
%   and computes the synthesis.
%
%   Note that U is of size [size(Y, 1), 3].

n = size(Y, 1);
dim = length(u);
assert(isvector(u));
assert(size(Y, 2) == dim);

U = zeros(n, 3);
parfor k=1:dim
    U = U + u(k) * squeeze(Y(:, k, :));
end

end