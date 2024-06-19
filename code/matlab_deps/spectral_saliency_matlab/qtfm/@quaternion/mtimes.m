function q = mtimes(l, r)
% *   Matrix multiply.
% (Quaternion overloading of standard Matlab function.)
 
% Copyright � 2005, 2009, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

ql = isa(l, 'quaternion');
qr = isa(r, 'quaternion');

if ql && qr
    
    % Both arguments are quaternions. There are four cases to handle,
    % dependent on whether the arguments are pure or full.
    
    pl = isempty(l.w);
    pr = isempty(r.w);
    
    if pl
        if pr
            % Both arguments are pure quaternions.
            q = l; % Create a result by copying one of the arguments. This
                   % avoids a costly call to the constructor.
            q.w = -(l.x * r.x + l.y * r.y + l.z * r.z);
            q.x =               l.y * r.z - l.z * r.y;
            q.y = - l.x * r.z             + l.z * r.x;
            q.z =   l.x * r.y - l.y * r.x;
        else
            % The right argument is full, but the left is pure.
            q = r;
            q.w = -(l.x * r.x + l.y * r.y + l.z * r.z);
            q.x =   l.x * r.w + l.y * r.z - l.z * r.y;
            q.y = - l.x * r.z + l.y * r.w + l.z * r.x;
            q.z =   l.x * r.y - l.y * r.x + l.z * r.w;
        end
    else
        if pr
            % The left argument is full, but the right is pure.
            q = l;
            q.w = - l.x * r.x - l.y * r.y - l.z * r.z;
            q.x =   l.w * r.x             + l.y * r.z - l.z * r.y;
            q.y =   l.w * r.y - l.x * r.z             + l.z * r.x;
            q.z =   l.w * r.z + l.x * r.y - l.y * r.x;
        else
            % Both arguments are full quaternions. The full monty.
            q = l;
            q.w =  l.w * r.w - (l.x * r.x + l.y * r.y + l.z * r.z);
            q.x =  l.w * r.x +  l.x * r.w + l.y * r.z - l.z * r.y;
            q.y =  l.w * r.y -  l.x * r.z + l.y * r.w + l.z * r.x;
            q.z =  l.w * r.z +  l.x * r.y - l.y * r.x + l.z * r.w;
        end
    end
   
else
    % One of the arguments is not a quaternion. If it is numeric, we can
    % handle it, and we must because the code above requires us to multiply
    % scalar parts by vector parts using a recursive call to this function.
    
    if ql && isa(r, 'numeric')
        q = l; % The left operand is a quaternion, so use it for the result
               % to avoid calling the constructor.     
        if ~isempty(l.w), q.w = q.w * r; end
        q.x = q.x * r; q.y = q.y * r; q.z = q.z * r;
    elseif isa(l, 'numeric') && qr
        q = r; % The right operand is a quaternion, so use it for the
               % result to avoid calling the constructor.  
        if ~isempty(r.w), q.w = l * q.w; end
        q.x = l * q.x; q.y = l * q.y; q.z = l * q.z;
    else
        error('Matrix multiplication of a quaternion by a non-numeric is not implemented.')
    end
end

% $Id: mtimes.m,v 1.9 2010/10/30 17:37:42 sangwine Exp $
