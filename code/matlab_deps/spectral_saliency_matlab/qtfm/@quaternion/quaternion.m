function q = quaternion(a0, a1, a2, a3)
% QUATERNION   Construct quaternions from components. Accepts the following
% possible arguments, which may be scalars, vectors or matrices:
%
% No argument     - returns an empty quaternion scalar, vector or matrix.
% One argument    - A quaternion argument returns the argument unmodified.
%                   A non-quaternion argument returns the argument in the
%                   scalar part and supplies a zero vector part with
%                   elements of the same type as the scalar part.
% Two arguments   - returns a quaternion, provided the first argument is
%                   numeric and the second is a pure quaternion.
% Three arguments - returns a pure quaternion scalar, vector or matrix,
%                   with an empty scalar part.
% Four arguments  - returns a full quaternion scalar, vector or matrix.

% Copyright � 2005, 2009, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargoutchk(0, 1, nargout)) % The number of input arguments is checked
                                 % as part of the switch statement below.
switch nargin
case 0

   % Construct an empty quaternion.

   q.w = []; q.x = []; q.y = []; q.z =[]; q = class(q, 'quaternion');

case 1

   if isa(a0, 'quaternion')
       q = a0; % a0 is already a quaternion, so return it.
   else
       % a0 is not a quaternion ...
       if isnumeric(a0)
           % ... but it is numeric. We can handle this, but we have to
           % supply the vector components of the quaternion, and they must
           % have the same type as a0 (e.g. double, int8), for which we use
           % the class option to the zeros function.
           
           t = zeros(size(a0), class(a0));
           q.w = a0; q.x = t; q.y = t; q.z = t; q = class(q, 'quaternion');
       else
           % ... or it isn't numeric.
           error(['Cannot construct a quaternion with a ',...
                  'non-numeric in the scalar part.']);
       end
   end

case 2
    
   % In this case, the first argument must be the scalar part, and thus
   % numeric, and the second must be the vector part, and thus a pure
   % quaternion (of the same component class as the first argument).

   if any(size(a0) ~= size(a1))
      error('Arguments must have the same dimensions')
   end

   if isnumeric(a0) && isa(a1, 'quaternion')
       c0 = class(a0); c1 = class(a1.x);
       if isempty(a1.w)
           if strcmp(c0, c1)
               % a1 is already a quaternion, so copy it to the output and
               % insert a0 as the w component.
               q = a1; q.w = a0;
           else
               error(['Classes of given scalar and vector parts', ...
                      ' differ: ', c0, ' ', c1])
           end
       else
           error('The second argument must be a pure quaternion.');
       end
   else
       error(['First argument must be numeric and the ',...
              'second must be a pure quaternion.']);
   end

case 3

   % Construct a pure quaternion, since we are given three arguments and
   % this is the only possibility. They must all be numeric and of the same
   % class.

   s0 = size(a0);
   if any(s0 ~= size(a1)) || any(s0 ~= size(a2))
       error('Arguments must have the same dimensions')
   end

   c0 = class(a0); c1 = class(a1); c2 = class(a2);
   % We test only the first argument for numeric, because if the classes
   % match, the other two must also be numeric.
   if isnumeric(a0) && strcmp(c0, c1) && strcmp(c1, c2)
       q.w = []; q.x = a0; q.y = a1; q.z = a2; q = class(q, 'quaternion');
   else
       error(['All three arguments must be numeric and of the', ...
              ' same class. Given: ', c0, ' ', c1, ' ', c2]);
   end

case 4 % Return a full quaternion. All four arguments must be numeric and
       % of the same class.

   s0 = size(a0);
   if any(s0 ~= size(a1)) || any(s0 ~= size(a2)) || any(s0 ~= size(a3))
       error('Arguments must have the same dimensions')
   end

   c0 = class(a0); c1 = class(a1); c2 = class(a2); c3 = class(a3);
   % We test only the first argument for numeric, because if the classes
   % match, the other three must also be numeric.
   if isnumeric(a0) && strcmp(c0, c1) && strcmp(c1, c2) && strcmp(c2, c3)
       q.w = a0; q.x = a1; q.y = a2; q.z = a3; q = class(q, 'quaternion');
   else
       error(['All four arguments must be numeric and of the', ...
              ' same class. Given: ', c0, ' ', c1, ' ', c2, ' ', c3]);
   end

otherwise
   error('Quaternion constructor takes from 0 to 4 arguments.');
end
end

% $Id: quaternion.m,v 1.14 2010/11/04 21:53:32 sangwine Exp $
