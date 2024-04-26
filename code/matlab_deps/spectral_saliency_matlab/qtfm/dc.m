function q = dc(A, B)
% DC: returns a quaternion from two complex numbers which are the
% Cayley-Dickson components of the quaternion. This function is the
% inverse of CD (q.v.).

% Copyright © 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(1, 2, nargin)), error(nargoutchk(0, 1, nargout))

if nargin == 1
    
    if ~isnumeric(A)
        error('Parameters must be numeric.')
    end

    if isa(A, 'quaternion')
        error('Parameters may not be quaternions.')
    end

    q = real(A) + imag(A) .* qi;
    
else
    
    % There must be two parameters.

    if any(size(A) ~= size(B))
        error('Parameters must be the same size.')
    end

    if ~isnumeric(A) || ~isnumeric(B)
        error('Parameters must be numeric.')
    end

    if isa(A, 'quaternion') || isa(B, 'quaternion')
        error('Parameters may not be quaternions.')
    end

    q = quaternion(real(A), imag(A), real(B), imag(B));
end

% $Id: dc.m,v 1.2 2009/02/08 19:18:16 sangwine Exp $

