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
function U = matrixU(dim, Z, Fc, V, ac)
%MATRIXU Creates the matrix U.
%
%   U = MATRIXU(dim, Z, Fc, V, ac) returns a symmetric full matrix U
%   with surface integrals u_{pq} = int_S Z_p*Z_q.
%
%   Note that size(U) = [dim, dim].
%
%   Note that U is symmetric and thus contains at most dim*(dim+1)/2
%   unique entries. However, using a sparse matrix representation with
%   adjusted vector-matrix multiplication results in slow matrix solves.
%   Thus, for the time being, a full matrix is used here to the
%   disadvantage of memory requirements.

U = zeros(dim, dim);
for p=1:dim
    for q=1:p
        U(p, q) = triangIntegral(Fc, V, Z(:, p) .* Z(:, q), ac);
        U(q, p) = U(p, q);
    end
end

end