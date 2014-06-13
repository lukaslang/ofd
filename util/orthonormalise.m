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
function [A, e, f] = orthonormalise(v, w)
%ORTHONORMALISE Computes an orthonormal basis.
%
%   [A, e, f] = orthonormalise(v, w) takes two matrices v, w of same size
%   representing a basis {v, w} of the tangent space and computes for each 
%   row a transformation matrix A such that {e, f} form an orthonormal 
%   basis and [e(k, :); f(k, :)] = A * [v(k, :); w(k, :)]. The basis is 
%   computed using the Gram-Schmidt process.
%
%   v and w are matrices of size n-by-3. Matrix A is of size 2-by-2-by-n
%   and e, f are of size n-by-3.
assert(size(v, 2) == 3);
assert(size(w, 2) == 3);
assert(size(v, 1) == size(w, 1));

% Number of points.
n = size(v, 1);

% Create coefficient matrix.
A = zeros(2, 2, n);

% Compute length of first basis vector.
len1 = sqrt(sum(v .^ 2, 2));

% Scale first basis vector to unit length.
e = bsxfun(@rdivide, v, len1);

% Compute inner product between v and w and subtract projection.
ip = sum(v .* w, 2);
f = w - bsxfun(@rdivide, bsxfun(@times, v, ip), len1 .^ 2);

% Compute length of second basis vector.
len2 = sqrt(sum(f .^ 2, 2));

% Scale second basis vector to unit length.
f = bsxfun(@rdivide, f, len2);

% Compute transformation matrix.
A(1, 1, :) = 1 ./ len1;
A(1, 2, :) = zeros(n, 1);
A(2, 1, :) = -ip ./ (len1 .^ 2) ./ len2;
A(2, 2, :) = 1 ./ len2;

end