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
function At = matrixAt(dim, Z, Fc, V, ac)
%MATRIXAT Creates the matrix \tilde{A}.
%
%   At = MATRIXAT(dim, Z, Fc, V, ac) returns a symmetric full matrix At
%   with surface integrals \tilde{a}_pq = int_S Z_p*Z_q.
%
%   Note that size(At) = [dim, dim].
%
%   Note that At is symmetric and thus contains at most dim*(dim+1)/2
%   unique entries. However, using a sparse matrix representation with
%   adjusted vector-matrix multiplication results in slow matrix solves.
%   Thus, for the time being, a full matrix is used here to the
%   disadvantage of memory requirements.

At = zeros(dim, dim);
for p=1:dim
    for q=1:p
        At(p, q) = triangIntegral(Fc, V, Z(:, p) .* Z(:, q), ac);
        At(q, p) = At(p, q);
    end
end

end