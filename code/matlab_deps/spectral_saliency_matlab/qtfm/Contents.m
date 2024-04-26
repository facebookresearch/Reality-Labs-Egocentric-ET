% Quaternion Toolbox (QTFM)
% Version 1.9 30-March-2011
% ---------------------------------------------------------------------
% Copyright (©) 2005-2011        Stephen J. Sangwine & Nicolas Le Bihan
% Email:        sjs@essex.ac.uk, nicolas.le-bihan@gipsa-lab.inpg.fr
%
% See the file Copyright.m for further details.
% ---------------------------------------------------------------------
%
% quaternion is a Matlab class library designed to extend Matlab in as
% natural a way as possible to handle quaternions. Many standard Matlab
% features have been implemented, including standard operators such as
% the arithmetic operators, matrix and elementwise products, indexing
% and indexed assignment using the colon operator, concatenation, end
% indexing, transpose and conjugate transpose, raising to a power, etc.
%
% A major objective has been to make it possible to write code that can
% work unchanged for real, complex, quaternion, and even complexified
% quaternion arrays, and for this reason standard Matlab functions have
% been overloaded for quaternion arrays.
%
% Most functions are fully vectorised and where they overload a standard
% Matlab function, the vectorisation is designed to work in the same way.
% Thus, for example, the power function (also accessible using the .^
% notation) will raise a single quaternion to a vector or matrix power in
% the same way as the Matlab power function. If the left and right
% arguments have the same size the power function operates elementwise,
% exactly as the standard Matlab function.
%
% A quaternion object, as implemented by this class, has a private
% implementation based on a structure array storing the four components
% of the quaternion. As with any standard Matlab type, a quaternion is by
% default a matrix. A single quaternion is simply a matrix of one element.
%
% Quaternions may be pure or full. Pure quaternions have no scalar part and
% operations that attempt to access the scalar part will result in an
% error, with the exception of the function scalar, which returns zero. The
% quaternions in a vector or matrix are either all pure or all full -- it
% is not possible to construct a vector or matrix with a mix of the two. A
% full quaternion with a scalar part which is zero differs from a pure
% quaternion, which has no scalar part.
%
% Quaternion matrices can be constructed either by using a constructor
% function, like this:
%
% q = quaternion(eye(5,5), zeros(5,5), randn(5,5), ones(5,5))
%
% or by using the three quaternion operators, named q1, q2 and q3, (these
% three operators are also available under the names qi, qj and qk) like this:
%
% q = eye(5,5) + zeros(5,5) * q1 + randn(5,5) * q2 + ones(5,5) * q3
%
% The components of the quaternion so constructed will have a type determined
% by the type supplied for the components (e.g. double, uint8, int16 ....).
% The components may be REAL or COMPLEX. A complex quaternion can also be
% constructed like this:
%
% q = complex(quaternion(1,2,3,4), 3)
%
% or using the Matlab value i:
%
% q = quaternion(1,2,3,4) + quaternion(5,6,7,8) .* i
%
% The components of a quaternion can be extracted using the functions:
%
%   scalar         - Scalar part of a quaternion.
%   vector, v      - Vector part of a quaternion (synonyms).
%   s, x, y, z     - Components of a quaternion.
%
%   real           - Real part of a (complex) quaternion.
%   imag           - Imaginary part of a (complex) quaternion.
%   complex        - Construct a complex quaternion from real quaternions.
%   quaternion     - Construct a quaternion from real or complex values.
%
% Quaternions may be decomposed into Cayley-Dickson form or constructed
% from Cayley-Dickson components using:
%
%   cd             - Decompose into two complex components.
%   dc             - Construct from two complex components.
%   cdpolar        - Decompose into two complex components, polar form.
%
% Other functions implemented are:
%
%   abs            - Modulus of a quaternion.
%   conj           - Conjugate (quaternion, complex or total).
%   unit           - Normalise a quaternion.
%   sign           - Equivalent to unit, cf Matlab sign.m.
%   inv            - Matrix and quaternion inverse.
%   axis           - Axis of a quaternion.
%   angle          - Angle or argument of a quaternion.
%
%   ceil, floor    - Round elements of quaternion towards plus or minus 
%   fix, round       infinity, zero, or nearest integer.
%
%   display        - Display array (does not show values)
%   disp           - Display without array name.
%   displayall     - Display components of quaternion.
%   show           - Shorter synonym for displayall.
%   char           - Convert quaternion to string.
%   fprintf        - Output quaternions to file.
%   write          - Write a quaternion array   to a text file.
%   read           - Read  a quaternion array from a text file.
%
%   cast, convert  - Convert components of quaternion to a different type.
%
%   scalar_product - Scalar (dot) product (elementwise).
%   cross          - Vector (cross) product (elementwise).
%   vector_product - Synonym for cross.
%
%   exp            - Exponential function.
%   log            - Natural logarithm.
%
%   sqrt           - Square root.
%
%   sin, cos, tan  - Trigonometric functions.
%   asin, acos,
%         atan     - Inverse trigonometric functions.
%   sinh, cosh,
%         tanh     - Hyperbolic functions.
%   asinh, acosh,
%          atanh   - Inverse hyperbolic functions.
%
%   blkdiag        - Block diagonal matrix.
%   diag           - Extract or construct a diagonal.
%   triu/tril      - Extract upper or lower triangular.
%   norm           - Vector and matrix norms.
%   sum            - Sum elements or columns.
%   prod           - Product of array elements.
%   cumsum         - Cumulative sum.
%   cumprod        - Cumulative product.
%   mean           - Mean of elements or columns.
%
%   eyeq           - Quaternion identity matrix.
%   onesq          - Quaternion matrix of ones.
%   zerosq, zerosv - Quaternion and pure quaternion matrices of zeros.
%   randv, randq   - Uniformly distributed unit vectors and quaternions.
%   randf          - Fisher distribution of unit vectors on sphere.
%   randvmf        - von Mises-Fisher distribution of unit quaternions.
%
%   ispure         - Test whether a quaternion (array) is pure.
%   isempty        - Test whether a quaternion (array) is empty.
%   isfinite       - Test whether a quaternion (array) is finite.
%   isinf          - Test whether a quaternion (array) is infinite.
%   isnan          - Test whether a quaternion (array) is NaN.
%   isreal         - Test whether a quaternion (array) is real.
%   ishermitian    - Test whether a quaternion (array) is Hermitian.
%   isunitary      - Test whether a quaternion (array) is unitary.
%
%   length         - Length of a quaternion vector.
%   size           - Size of a quaternion array.
%   ndims          - Number of array dimensions.
%   numel          - Number of elements in a quaternion array. NB This does
%                    not return prod(size(q)) as might be expected.
%   repmat         - Replicate and tile a quaternion array.
%   reshape        - Change size of array.
%   squeeze        - Remove singleton array dimensions.
%   cat            - Concatenate arrays.
%   permute        - Permute dimensions of N-D array.
%   ipermute       - Inverse permute dimensions of N-D array.
%
%   det            - Determinant.
%   svd            - Singular value decomposition.
%   eig            - Eigenvalue decomposition.
%   lu             - LU decomposition.
%   qr             - QR decomposition.
%
%   expm           - Matrix exponential.
%   logm           - Matrix logarithm.
%   sqrtm          - Matrix square root.
%
%   slerp          - Spherical linear interpolation.
%
%   imreadq        - Read RGB image into quaternion matrix.
%   imwrite        - Write quaternion matrix to RGB image.
%   image          - Display quaternion array as an image.
%   scatter3       - 3D scatter plot of pure quaternion vector.
%
%   adjoint        - The adjoint matrix of a quaternion matrix.
%   unadjoint      - Construct a quaternion matrix from its adjoint.
%
%   conv, conv2    - Convolution.
%
%   fft            - One dimensional (default) quaternion Fourier transform.
%   fft2           - Two dimensional (default) quaternion Fourier transform.
%   qfft           - One dimensional left or right one-dimensional QFFT.
%   qdft           - One dimensional left or right one-dimensional QDFT.
%   qfft2          - Two dimensional left or right two-dimensional QFFT.
%   qdft2          - Two dimensional left or right two-dimensional QDFT.
%   fftshift       - Quaternion overloading of the standard Matlab function.
%   .............  - All of the above have inverses, prefixed with 'i'.
%
%   bsxfun         - Binary Singleton Expansion Function.
%
% The following builtin Matlab functions also work for quaternion arrays:
%
%   flipud, fliplr, flipdim, rot90, circshift, trace
%   isequal, isscalar, isvector
%
% The following Matlab functions also work for quaternion arrays:
%
%   cov, dot, rank, std, var (there may be others)
%
% There are some auxiliary functions which are used to compute some of the
% more elaborate functions above, such as svd, eig and qdft. These are:
%
% change_basis, orthonormal_basis, orthogonal.
% householder_vector, householder_matrix, bidiagonalize, tridiagonalize.
%
% Unimplemented functions (for which a placeholder file is provided) are:
%
% arrayfun       - Apply a function to each element of an array.
% funm           - Matrix function.
%
% These unimplemented functions may be implemented at a later date.
%
% Some test code is provided in the directory 'test'. To run it set the
% working directory to 'test' and type 'test'. This runs all the test
% code. Alternatively, there is a menu item accessible from the Matlab
% Start menu which will run the test code.
%
% For more information, use help, accessible from the command window,
% the Matlab Start menu (Toolboxes -> Quaternion), or the Matlab help
% browser (Quaternion Toolbox).

% $Id: Contents.m,v 1.43 2011/01/09 18:09:30 sangwine Exp $
