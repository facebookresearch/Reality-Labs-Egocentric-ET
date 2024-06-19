function r = eq(a, b)
% ==  Equal.
% (Quaternion overloading of standard Matlab function.)
%
% If one of the operands is not a quaternion and the other has zero vector part,
% the result is obtained by comparing the non-quaternion operand with the scalar
% part of the quaternion operand.

% Copyright � 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout)) 

if isa(a, 'quaternion') && isa(b, 'quaternion')
    
    pa = isempty(a.w);
    pb = isempty(b.w);
    if pa && pb
        r = a.x == b.x & a.y == b.y & a.z == b.z;
    elseif pa
        % a is pure, but b isn't, so they can only be equal if the
        % scalar part of b is zero and the vector parts are equal.
        
        r = b.w == 0 & vector(a) == vector(b);
    elseif pb
        % b is pure, but a isn't, so they can only be equal if the
        % scalar part of a is zero and the vector parts are equal.
        
        r = a.w == 0   & vector(a) == vector(b);
    else
        r = a.w == b.w & vector(a) == vector(b);
    end
    
else
    % One of the arguments is not a quaternion (the other must be,
    % otherwise Matlab would not call this function).
    
    if isa(a, 'quaternion')
        if isempty(a.w)
            % a has no scalar part, and b has no vector part, since it is not
            % a quaternion. The result can only be true if b is zero, and a
            % has zero vector part.
            
            r = b == 0 & a == quaternion(0, 0, 0);
        else
            % a is a full quaternion, so the result can be true only if a has
            % zero vector part, and b is equal to the scalar part of a.
            
            r = a.w == b & vector(a) == quaternion(0, 0, 0);
        end
    elseif isa(b, 'quaternion')
        r = b == a; % Swap the order and compare them using the code above.
    else
        error('Programming error: neither argument is a quaternion?');    
    end
end

% $Id: eq.m,v 1.7 2010/11/04 21:50:45 sangwine Exp $
