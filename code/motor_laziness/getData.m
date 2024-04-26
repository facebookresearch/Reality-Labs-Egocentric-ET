function [d3] = getData(directory)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% getData.m

d = load('data.mat');
d = d.TIMESERIES;

for iter = 1:2 % iterate to output two datasets, one for GazeInHead, the second for GazeInSpace
    
    if iter == 1
        gazeData = {d.GazeInHeadX d.GazeInHeadY};
    elseif iter == 2
        gazeData = {d.GazeInSpaceX d.GazeInSpaceY};
    end
    
    
    %%%%%%%%%%%% downsample IMU data to sampling rate of gaze data
    for kk = 1:length(d.Class) % loop through recordings
        
        d.HeadVelAngMag2{kk} = resample(d.HeadVelAngMag{kk},1,2); % downsampling from 100 to 50 Hz
        
        % then trimming excess or padding with nans
        if length(gazeData{1}{kk}) > length(d.HeadVelAngMag2{kk})
            temp123 = NaN(length(gazeData{1}{kk}),1);
            temp123(1:length(d.HeadVelAngMag2{kk})) = d.HeadVelAngMag2{kk};
            d.HeadVelAngMag2{kk} = temp123;
        else
            d.HeadVelAngMag2{kk} = d.HeadVelAngMag2{kk}(1:length(gazeData{1}{kk}));
        end
        
    end
    
    %%%%%%%%%%%%%% estimate fixation locations and saccades
    for kk = 1:length(d.Class) % loop through recordings
        fixationOrNot = d.Class{1,kk};
        v1 = find(fixationOrNot==1);
        v2 = diff(v1);
        indsOfLastOne = find(v2>1); % ind in fixation segment where there is the final 1.
        clear indsOfFixations; clear fixLocation; clear fixDur; clear headSpeed;
        for ii = 1:length(indsOfLastOne) %% ii is number of fixations
            if ii == 1
                indsOfFixations{ii} = v1(1:indsOfLastOne(ii));
            else
                indsOfFixations{ii} = v1(indsOfLastOne(ii-1)+1:indsOfLastOne(ii));
            end
            
            % Now avg over gaze coordinate during detected fixation to get fixation
            % location in deg. This is the scanpath:
            fixLocation(ii,:) = nanmean([gazeData{1}{1, kk}(indsOfFixations{ii}) gazeData{2}{1, kk}(indsOfFixations{ii})]);
            saccsCart = diff(fixLocation);
            
            fixDur(ii,:) = length(gazeData{1}{1, kk}(indsOfFixations{ii}))*.02; % fixation duration in seconds
            
            
            % get median head speed during these fixation intervals as well
            headSpeed(ii,:) = median(d.HeadVelAngMag2{kk}(indsOfFixations{ii}));
            
        end
        
        d2.fixLocationX{kk} = fixLocation(:,1); d2.fixLocationY{kk} = fixLocation(:,2);
        d2.saccsCart{kk} = saccsCart;
        d2.fixDur{kk} = fixDur;
        
        d2.distFromCenterOfFOV{kk} = vecnorm([d2.fixLocationX{kk} d2.fixLocationY{kk}]')';
        
        d2.headSpeed{kk} = headSpeed;
        
    end
    
    %%%%%%%%%%% get metadata about task and subject for each recording
    
    fff = cat(1,d.Task{:});
    d2.taskList = unique(fff);
    
    for bb = 1:length(d2.taskList)
        d2.indicesForTaskbb{bb} = find(strcmp(fff, d2.taskList{bb}));
        d2.subListForTaskbb{bb} = d.SubID(find(strcmp(fff, d2.taskList{bb})));
    end
    
    
    %%%%%%%%%%% Extract saccades for each task
    d2.saccsPooledByTask = cell(length(d2.taskList),1);
    for bb = 1:length(d2.taskList)    
        for jj = 1:length(d2.indicesForTaskbb{bb})
            curInd = d2.indicesForTaskbb{bb}(jj);
            d2.saccsPooledByTask{bb} = [d2.saccsPooledByTask{bb}; d2.saccsCart{curInd}];
        end
    end
    
    %%%%%%%%%%%% save out raw gaze for calculating gaze dispersion for Fig. 4
    if iter == 1
        d2.GazeX = d.GazeInHeadX;
        d2.GazeY = d.GazeInHeadY;
    elseif iter ==2
        d2.GazeX = d.GazeInSpaceX;
        d2.GazeY = d.GazeInSpaceY;
    end
    
    %%%%%%%%%% now output gazeData
    if iter == 1
        d3.GazeInHead = d2;
    elseif iter == 2
        d3.GazeInSpace = d2;
    end
    
    
end


end