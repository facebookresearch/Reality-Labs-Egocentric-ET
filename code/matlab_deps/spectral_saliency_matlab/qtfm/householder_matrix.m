function [h] = householder_matrix(a, v)
% h = householder_matrix(a, v) returns the Householder matrix that
% will zero all elements of a except those corresponding to (any)
% non-zero elements of v. If v has zero norm, an identity matrix is
% returned. The matrix returned is either a left or right Householder
% matrix dependent on whether the input parameters are column or row
% vectors respectively.

% Copyright � 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Reference:
%
% Sangwine, S. J. and Le Bihan, N.,
% Quaternion singular value decomposition based on bidiagonalization
% to a real or complex matrix using quaternion Householder transformations,
% Applied Mathematics and Computation, 182(1), 1 November 2006, 727-738, 
% DOI:10.1016/j.amc.2006.04.032.
%
% S. J. Sangwine and N. Le Bihan,
% Quaternion Singular Value Decomposition based on Bidiagonalization
% to a Real Matrix using Quaternion Householder Transformations,
% arXiv:math.NA/0603251, 10 March 2006. Available at http://www.arxiv.org/.

error(nargchk(2, 2, nargin)), error(nargoutchk(1, 1, nargout))

if size(a) ~= size(v)
    error('Input parameters must be vectors of the same size, both row, or both column')
end

[m, n] = size(a);

col_vector = n == 1 & m == length(a);
row_vector = m == 1 & n == length(a);

if col_vector == 0 && row_vector == 0
    error('Parameters must be either column or row vectors.');    
end

h = quaternion(eye(length(a)));

n = norm(v); % If v has zero norm, we return an identity matrix.

if n ~= 0
    [u, zeta] = householder_vector(a, v ./ n);
    if size(u) ~= size(a)
        error('Result from function householder_vector is not of correct type (row, column)');
    end
    if col_vector
        h = (1 ./ zeta) .* (h - u * u');
    else
        h = (h - u.' * conj(u)) .* (1 ./ zeta);
    end
end

% $Id: householder_matrix.m,v 1.4 2009/02/08 19:18:16 sangwine Exp $

