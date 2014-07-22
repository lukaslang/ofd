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
function test_suite = lebedevIntegralTest
    initTestSuite;
end

function constantFunctionTest

% Create constant function.
f = @(x, y, z) x.^2 + y.^2 + z.^2;

% Compute surface integral over function which is constant one.
v = lebedevIntegral(f, 3);

assertTrue(isscalar(v));
assertAlmostEqual(v, 4*pi, 1e-10);

end

function constantZeroFunctionTest

% Create constant function.
f = @(x, y, z) zeros(length(x), 1);

% Compute surface integral over function which is constant one.
v = lebedevIntegral(f, 3);

assertTrue(isscalar(v));
assertAlmostEqual(v, 0, 1e-10);

end

function constantSpharmTest

% Create spherical harmonic of degree zero.
f = @(x, y, z) spharmn(0, [x, y, z]).^2;

% Compute scalar product of spherical harmonic.
v = lebedevIntegral(f, 3);

assertTrue(isscalar(v));
assertAlmostEqual(v, 1, 1e-10);

end

function spharmOnbOrderTest

for N=1:5
    % Create spherical harmonics of degree N.
    f = @(x, y, z) spharmn(N, [x, y, z]);
    % Extract order k.
    Yk = @(f, k) f(:, k);
    
    % Run through all orders and test orthonormality with respect to L2 
    % inner product.
    for k=1:2*N+1
        for l=1:2*N+1
            % Compute scalar product of spherical harmonic.
            d = @(x, y, z) Yk(f(x, y, z), k) .* Yk(f(x, y, z), l);
            v = lebedevIntegral(d, 11);
            assertTrue(isscalar(v));
            assertAlmostEqual(v, double(k == l), 1e-10);
        end
    end
end
end

function spharmOnbDegreeTest

M = 5;
% Create spherical harmonics of degree M.
fm = @(x, y, z) spharmn(M, [x, y, z]);

N = 4;
% Create spherical harmonics of degree N.
fn = @(x, y, z) spharmn(N, [x, y, z]);

% Extract order k.
Yk = @(f, k) f(:, k);

% Run through all orders and test orthonormality with respect to L2 inner
% product.
for k=1:2*M+1
    for l=1:2*N+1
        % Compute scalar product of spherical harmonic.
        d = @(x, y, z) Yk(fm(x, y, z), k) .* Yk(fn(x, y, z), l);
        v = lebedevIntegral(d, 11);
        assertTrue(isscalar(v));
        assertAlmostEqual(v, 0, 1e-10);
    end
end
end