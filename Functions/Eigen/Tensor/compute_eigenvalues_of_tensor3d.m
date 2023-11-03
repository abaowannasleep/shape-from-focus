function [lambda1 lambda2 lambda3] = compute_eigenvalues_of_tensor3d(t11, t12, t13, t22, t23, t33)
% COMPUTE_EIGENVALUES_OF_TENSOR3D Estimate the eigenvalues of the real symmetric 3x3 matrix T
%
% [lambda1 lambda2 lambda3] = compute_eigenvalues_of_tensor3d(t11, t12, t13, t22, t23, t33)
%

% Copyright (c) 2012 Daniel Forsberg
% danne.forsberg@outlook.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Taken from http://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
% Matrix A is replaced with the elements t11, t12, t13, t22, t23 and t33.
% Also added part to always make sure that eigenvalues a positive

q = (t11 + t22 + t33)/3;
p = (t11 - q).^2 + (t22 - q).^2 + (t33 - q).^2 + 2 * (t12.^2 + t13.^2 + t23.^2);
p = sqrt(p / 6) + eps;
B11 = (1 ./ p) .* (t11 - q);
B12 = (1 ./ p) .* (t12);
B13 = (1 ./ p) .* (t13);
B22 = (1 ./ p) .* (t22 - q);
B23 = (1 ./ p) .* (t23);
B33 = (1 ./ p) .* (t33 - q);

r = (- B33.*B12.^2 + 2*B12.*B13.*B23 - B22.*B13.^2 - B11.*B23.^2 + B11.*B22.*B33) / 2;

% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
phi = acos(r) / 3;
phi(r <= -1) = pi/3;
phi(r >= 1) = 0;

% the eigenvalues satisfy lambda3 <= lambda2 <= lambda1
lambda1 = q + 2 * p .* cos(phi);
lambda3 = q + 2 * p .* cos(phi + pi * (2/3));
lambda2 = 3 * q - lambda1 - lambda3;     % since trace(A) = lambda1 + lambda2 + lambda3

ind = find(lambda3 < 0);
if ~isempty(ind)
    lambda1(ind) = lambda1(ind) - 2*lambda3(ind);
    lambda2(ind) = lambda2(ind) - 2*lambda3(ind);
    lambda3(ind) = lambda3(ind) - 2*lambda3(ind);
end
 No newline at end of file
