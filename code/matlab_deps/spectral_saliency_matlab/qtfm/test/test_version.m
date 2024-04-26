function test_version
% Test code to verify the Matlab or GNU Octave version.

% Copyright © 2008, 2010, 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Checking Matlab version ...');

V = ver('Matlab'); % This gets a structure containing information about the
                   % Matlab version, if we are running under Matlab. If we
                   % are not, the result will be empty.
if ~isempty(V)

    S = V.Version; % This gets and displays a string containing the Matlab
    disp(S);       % version. S starts with 7.7 or 7.10 or similar.
    
    M = 'Matlab version must be 7.4 or later';
    N = 'Cannot make sense of Matlab version string';

    if str2double(S(1)) < 7
        error(M)
    end
    
    % TODO Clean up the parsing of the version string. It must be possible
    % to extract the major and minor numbers by looking for the dots,
    % rather than messing with the positions of the digits.
    L = length(S);
    if L == 3
        % The version number is of the form '7.x'.
        R = str2double(S(3));
    elseif L == 4
        % The version number is of the form '7.xy'.
        R = str2double(S(3:4));
    elseif L == 6
        % The version number is of the form '7.xy.z' (first used with
        % R2010bSP1).
        R = str2double(S(3:4));
    else
        error(N);
    end
    if R < 4
        % We need 7.4 minimum because of the assert function and also the
        % bsxfun function needed to support var, std and cov.
        error(M) % We must have version 7.3 or earlier.
    end
    
    disp('Matlab version is OK.');

else
    
    disp('Not running under Matlab, checking for Octave ...');
    
    V = ver('Octave');
    if isempty(V)
        error('Not running under Matlab or Octave, giving up!');
    end
    
    S = V.Version;
    if S(2) == '.' && S(4) == '.'
        if str2double(S(1)) < 3
           error('Minimum major version required is 3.');
        else
            if str2double(S(1)) == 3 && str2double(S(3)) < 2
                error('Minimum minor version required is 2.');
            end
        end
    else
        error('Cannot make sense of Octave version string.')
    end

    disp('Octave version is OK.');
    
end

% $Id: test_version.m,v 1.9 2011/03/25 16:23:14 sangwine Exp $
