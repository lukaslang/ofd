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
function test_suite = integrateflow2Test
    initTestSuite;
end

function resultTest   
    h = 0.1;
    maxit = 1;
    P = [1, 1;
        0, 1;
        1, 0;
        0, 0];
    v = ones(4, 2);
    s = [0, 0];
    verts = integrateflow2(P, v, s, h, maxit);
    assertEqual(length(verts), 1);
    assertEqual(verts{1}, [0, 0; 0.1, 0.1]);
end