function Q = imreadq(varargin)
% IMREADQ  Read image from graphics file and return quaternion image array.
%
% This function takes the same parameters as the Matlab function IMREAD,
% but it returns a quaternion array, with elements of type uint8 or uint16
% depending on the pixel data in the image file (8-bit or 16-bit). The way
% the image pixel data is converted to quaternions depends on the type of
% image stored in the file as follows:
%
% 1. RGB colour image in file. The data is returned as a pure quaternion
%    array with the RGB components in the XYZ components of the quaternion.
%
% No other cases are handled at present. For more precise control over the
% conversion from image to quaternion array, use the Matlab function IMREAD
% and then code the conversion explicitly.

% Copyright � 2009 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(nargoutchk(0, 1, nargout)) % The number of input parameters is
                                 % checked by imread, so we don't check it
                                 % here.
                                 
A = imread(varargin{:}); % Read the specified image (call to standard
                         % Matlab function).
                      
% Now convert the data in A into quaternion format depending on what it is.

N = size(A, 3); % Find out how many image components are in the array A.

if N == 3 % The image is RGB.
    Q = quaternion(A(:,:,1), A(:,:,2), A(:,:,3));
else
    error('The number of components in the image is not 3. Not handled.');
end

% $Id: imreadq.m,v 1.1 2009/03/02 15:12:02 sangwine Exp $
