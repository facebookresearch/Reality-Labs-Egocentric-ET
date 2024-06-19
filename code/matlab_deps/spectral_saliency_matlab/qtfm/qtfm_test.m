function qtfm_test
% Run qtfm test code.
%
% Copyright � 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

current_dir = pwd;

root = [qtfm_root filesep 'test'];
cd(root)
test

cd(current_dir);

% This file is provided so that test code can be run from the Start menu.
% It would be more elegant if the test code could be called directly, but
% the problem is how to specify the location in the info.xml file, since it
% will vary according to where the toolbox has been installed. Provided the
% path is set up correctly, meaning that this file can be found on the
% path, it will run the test code, because the cd command operates relative
% to the current directory, and the qtfm_root function determines where the
% QTFM files have been installed (all is that is necessary is that
% qtfm_root.m is on the current path).
% $Id: qtfm_test.m,v 1.2 2009/02/08 19:18:16 sangwine Exp $

