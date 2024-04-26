function [e] = fig2(d) 
% Copyright (c) Meta Platforms, Inc. and affiliates.
% fig2
% save out samples from model simulations for other figures

d = d.GazeInHead;

%%  DATA: bootstrapping of saccade vectors to equalize N across tasks ****** %%

clear saccsPooledByTaskBootstrap

% find the N from the task with the highest N, that will be the number of
% fixations we draw from each task to equalize data across tasks
for ii = 1:length(d.taskList); xvx(ii) = size(d.saccsPooledByTask{ii},1); end
maxN = max(xvx);

for ii = 1:length(d.taskList)
    [theta1,rho1] = cart2pol(d.saccsPooledByTask{ii}(:,1),d.saccsPooledByTask{ii}(:,2)); % in this coordinate system, up is +90, right is 0, left, is +-180, down is -90
    theta1 = -theta1;
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
    both = [relativeSaccAngleDiff1 relativeSaccAmpDiff1];
    
    bs_inds = randi(size(both,1),maxN,1);
    saccsPooledByTaskBootstrap{ii} = both(bs_inds,:);
end

allSaccDiffs = cat(1,saccsPooledByTaskBootstrap{:});

relativeSaccAngleDiff1 = allSaccDiffs(:,1);
relativeSaccAmpDiff1 = allSaccDiffs(:,2);

%% Plot pooled tasks relative saccade joint distribution

% need to now wrap these around, so angles greater than 180º or less than
% -180º end up in right place
relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);

figure(2); subplot(1,2,1);
histogram2Polar(relativeSaccAngleDiff1,relativeSaccAmpDiff1+6,.5); % this doesn't work because the plot won't display negative radii
vv = [relativeSaccAngleDiff1,relativeSaccAmpDiff1];
title('Data')
colorbar

ratioReturnToForewardSaccAllSaccs = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );

disp(['Data: Ratio of Prevalence of Return To Forward Saccades = ',num2str(ratioReturnToForewardSaccAllSaccs)]);





%% MODEL: Now do the same but for the model


%% First compute fixation probability and duration maps for data, for each task separately

xBins = linspace(-40,40,60);
yBins = linspace(-25,25,60);

for bb = 1:length(d.indicesForTaskbb)
    
    fixPosPerTaskX = cat(1,d.fixLocationX{d.indicesForTaskbb{bb}});
    fixPosPerTaskY = cat(1,d.fixLocationY{d.indicesForTaskbb{bb}});
    fixDurPerTask = cat(1,d.fixDur{d.indicesForTaskbb{bb}});
    
    for ii = 1:length(xBins)-1
        for jj = 1:length(yBins)-1
            clear justDataInBin;
            justDataInBin = fixDurPerTask(fixPosPerTaskX>=xBins(ii) & fixPosPerTaskX<xBins(ii+1) & fixPosPerTaskY>=yBins(jj) & fixPosPerTaskY<yBins(jj+1));
            meanmapPerTask{bb}(jj,ii) = nanmean(justDataInBin);
            nPointsmapPerTask{bb}(jj,ii) = length(justDataInBin);
            
            % uncomment if you want to threshold so that bins with below 5
            % datapoints are removed. this way the colormap limits will be
            % more accurate
            if nPointsmapPerTask{bb}(jj,ii) < 5
                meanmapPerTask{bb}(jj,ii) = NaN;
            end
            
        end
    end
end

%% Now sample fixations from those distributions

colorsss = hsv(9);
xBins = linspace(-40,40,59);
yBins = linspace(-25,25,59);

% find the N from the task with the highest N, that will be the number of
% fixations we draw from each task to equalize data across tasks
for ii = 1:length(d.taskList); xvx(ii) = size(d.saccsPooledByTask{ii},1); end
maxN = max(xvx);

for bb = 1:length(d.indicesForTaskbb)

    totalNPoints = maxN;
    pMap = nPointsmapPerTask{bb}./totalNPoints;

    for ii = 1:totalNPoints

        [temp1, temp2] = pinky(xBins,yBins,pMap);
        e.locSampleT{bb}(ii,:) = [temp1 temp2];

        closestIndex1 = find(xBins==temp1);
        closestIndex2 =  find(yBins==temp2);

        e.durSampleT{bb}(ii) =  meanmapPerTask{bb}(closestIndex2,closestIndex1);
        e.distFromFOVcenterSampleT{bb}(ii) = vecnorm(e.locSampleT{bb}(ii,:)');

    end 
end

% Note: sometimes the duration sample will be a NaN. this is normal,
% because the fixation duration maps have no data in certain locations by
% the virtue of the fact that no one ever look there in the experimental
% data. It's better to treat this as missing data in the simulations, and
% just average over the simulated data that we have actual information
% about.


% now compute saccades and pool over tasks to plot model simulation of relative saccade joint distribution

    allModelPos = cat(1,e.locSampleT{:});
    saccsCartSample = diff(allModelPos);

    [theta1,rho1] = cart2pol(saccsCartSample(:,1),saccsCartSample(:,2));
    theta1 = -1*theta1;
    
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as Bays and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);
    
    % remove  inf values from the log transform
    relativeSaccAmpDiff1noInf = relativeSaccAmpDiff1;
    relativeSaccAngleDiff1noInf = relativeSaccAngleDiff1;
    infInds = (relativeSaccAmpDiff1noInf==inf | relativeSaccAmpDiff1noInf==-inf);
    relativeSaccAmpDiff1noInf(infInds) = [];
    relativeSaccAngleDiff1noInf(infInds) = [];
    
    figure(2); subplot(1,2,2)
    histogram2Polar(relativeSaccAngleDiff1noInf,relativeSaccAmpDiff1noInf+6,.5, 'RTicks',[]); % this doesn't work because the plot won't display negative radii
    title('Model');
    colorbar
    
    vv = [relativeSaccAngleDiff1noInf,relativeSaccAmpDiff1noInf];
    
    ratioReturnToForewardSaccAllSaccs = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );
    disp(['Model: Ratio of Prevalence of Return To Forward Saccades = ',num2str(ratioReturnToForewardSaccAllSaccs)]);

   
end
