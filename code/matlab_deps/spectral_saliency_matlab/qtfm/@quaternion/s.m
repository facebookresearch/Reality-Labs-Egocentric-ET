function x = s(q)
% S(Q) Scalar part of a full quaternion.

% Copyright © 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% error(nargchk(1, 1, nargin)), error(nargoutchk(0, 1, nargout))

x = q.w;

% $Id: s.m,v 1.3 2010/05/13 20:34:25 sangwine Exp $
