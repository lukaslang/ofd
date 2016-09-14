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
function p = projecttoplane(v)
%PROJECTTOPLANE Projects a vector field to the plane and rescales it to the
%original length.
%
%   p = PROJECTTOPLANE(v) takes a vector field v in R^3 and projects it to
%   R^2. In addition, it is rescaled to the original length.
%
%   v must be of size [n, 3], where n is the number of vectors. p is of
%   equal size.

assert(size(v, 2) == 3);

% Project to R^2.
p = v;
p(:, 3) = 0;

% Compute vector lengths.
lenv = sqrt(sum(v .^2, 2));
lenp = sqrt(sum(p .^2, 2));

% Find indices of vectors to be rescaled.
idx = lenv > 0 & lenp > 0;
p(~idx, :) = 0;

% Rescale projected vector.
p(idx, :) = bsxfun(@times, p(idx, :), lenv(idx, :) ./ lenp(idx, :));

end