function test_fft
% Test code for the fast quaternion Fourier transform.

% This code tests the following functions:
%
%  fft  fft2
% ifft ifft2

% Copyright � 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% We have to test the 1D and 2D functions separately because they do not
% call each other. In addition, we need to verify the results of the FFTs
% against the DFT code in the corresponding qdft functions. The fft and
% fft2 functions use a default axis, whereas the qdft and qdft2 functions
% don't, so it is important below that we call the dft functions with the
% same axis as the fft defaults (see the private function dft_axis). There
% seems to be no way to call dft_axis from here, even though we know the
% relative path or can construct an absolute path.

% Test 1. Real data, and therefore real axis, 1D and 2D.
% Check that the transform inverts correctly and then that it agrees with
% the DFT code.

disp('Testing fast Fourier transform functions ...')

T = 1e-12;

q = randq(10) .* randn(10);

RA = unit(quaternion(1,1,1)); % Real axis.

compare(q,  ifft(fft(q)),            T, 'fft/ifft failed test 1.');
compare(q, iqdft(fft(q), RA, 'L'),   T, 'fft/iqdft failed test 1.');
compare(q, ifft2(fft2(q)),           T, 'fft2/ifft2 failed test 1.');
compare(q, iqdft2(fft2(q), RA, 'L'), T, 'fft2/iqfft2 failed test 1.');

% Test 2. Complex data, and therefore complex axis, 1D and 2D.
% Check that the transform inverts correctly and then that it agrees with
% the DFT code.

b = complex(randq(10), randq(10)) .* randn(10);

CA = quaternion(1,1,1) + quaternion(0,1,-1).*1i; % Complex axis.

compare(b, ifft(fft(b)),             T, 'fft/ifft failed test 2.');
compare(b, iqdft(fft(b), CA, 'L'),   T, 'fft/iqdft failed test 2.');
compare(b, ifft2(fft2(b)),           T, 'fft2/ifft2 failed test 2.');
compare(b, iqdft2(fft2(b), CA, 'L'), T, 'fft2/iqfft2 failed test 2.');

disp('Passed');

% $Id: test_fft.m,v 1.5 2010/05/12 21:08:12 sangwine Exp $
