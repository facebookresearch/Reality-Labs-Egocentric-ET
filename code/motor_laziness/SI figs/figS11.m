function [] = figS11(d)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS11

d = d.GazeInHead;

%% DATA: boostrap to get equal data across tasks

% only preserve samples where headspeed is below 3 deg/s
for kk = 1:length(d.fixDur)
    fixPosX2{kk} = d.fixLocationX{kk}(d.headSpeed{kk}<3);
    fixPosY2{kk} =  d.fixLocationY{kk}(d.headSpeed{kk}<3);
    fixDur2{kk} =  d.fixDur{kk}(d.headSpeed{kk}<3);
    distFromCenterOfFOV2{kk} = d.distFromCenterOfFOV{kk}(d.headSpeed{kk}<3);
    
    gazeX2{kk} = d.GazeX{kk}(d.headSpeed{kk}<3);
    gazeY2{kk} =  d.GazeY{kk}(d.headSpeed{kk}<3);
end


for bb = 1:length(d.indicesForTaskbb)  
    fixPosPerTaskX{bb} = cat(1,fixPosX2{d.indicesForTaskbb{bb}});
    fixPosPerTaskY{bb} = cat(1,fixPosY2{d.indicesForTaskbb{bb}});
    fixDurPerTask{bb} = cat(1,fixDur2{d.indicesForTaskbb{bb}});
    distFromCenterOfFOVPerTask{bb} = cat(1,distFromCenterOfFOV2{d.indicesForTaskbb{bb}});
    
    headSpeedPerTask{bb} = cat(1,d.headSpeed{d.indicesForTaskbb{bb}});
end


% now bootstrap!

% find the N from the task with the highest N, that will be the number of
% fixations we draw from each task to equalize data across tasks
for ii = 1:length(d.taskList); xvx(ii) = size(fixPosPerTaskX{ii},1); end
maxN = max(xvx);

for ii = 1:length(d.taskList)
    
    both = [fixPosPerTaskX{ii} fixPosPerTaskY{ii}];
    
    bs_inds = randi(size(both,1),maxN,1);
    fixPosPerTaskBootstrap{ii} = both(bs_inds,:);
    
    fixDurPerTaskBootstrap{ii} = fixDurPerTask{ii}(bs_inds);
    
end
% remove walk task because the amount of data pre-bootstrapping is way too
% low,- 16 datapoints, which will just overwhelm pooled data with noise.
% ***
fixDurPerTaskBootstrap{9}=[];
fixPosPerTaskBootstrap{9} = [];

% pool data across tasks
posAlllll2 = cat(1,fixPosPerTaskBootstrap{:});
fixAllll2= cat(1,fixDurPerTaskBootstrap{:});

% spatial binning
xBins = linspace(-40,40,60);
yBins = linspace(-25,25,60);

% generate fixation position and duration maps
for ii = 1:length(xBins)-1
    for jj = 1:length(yBins)-1
        justDataInBin = fixAllll2(posAlllll2(:,1)>=xBins(ii) & posAlllll2(:,1)<xBins(ii+1) & posAlllll2(:,2)>=yBins(jj) & posAlllll2(:,2)<yBins(jj+1));
        meanmap1(jj,ii) = nanmean(justDataInBin);
        semMap(jj,ii) = nanstd(justDataInBin)./sqrt(length(justDataInBin));
        nPointsmap1(jj,ii) = length(justDataInBin);
    end
end


figure;
subplot(3,2,1);
imagesc(xBins,yBins,nPointsmap1)
axis equal; hline(0,'-k'); vline(0,'-k');
set(gca,'YDir','normal')
colormap(flipud(parula))
xlim([-40 40])
ylim([-25 25])
colorbar
xlabel('Horizontal eye-in-head position (º)')
ylabel('Vertical gaze position (º)')

subplot(3,2,2);
imagesc(xBins,yBins,meanmap1,'AlphaData',nPointsmap1>=11 ) % 11 obeys the ~ 1.12 data:bootstrapped&headMvmtThresholdedData ratio, equivalent to removing 10 datapoints
axis equal; hline(0,'-k'); vline(0,'-k');
set(gca,'YDir','normal')
colormap(flipud(parula))
caxis([.1 .8])
ylim([-25 25])
xlim([-40 40])
colorbar
xlabel('Horizontal eye-in-head position (º)')
ylabel('Vertical gaze position (º)')



%% Now do fig 2 and 3 equivalent but only for model, 


%% Fig 2 equiv. 
xBins = linspace(-40,40,60);
yBins = linspace(-25,25,60);

for bb = 1:length(d.indicesForTaskbb)
   
    for ii = 1:length(xBins)-1
        for jj = 1:length(yBins)-1
            clear justDataInBin;
            justDataInBin = fixDurPerTask{bb}(fixPosPerTaskX{bb}>=xBins(ii) & fixPosPerTaskX{bb}<xBins(ii+1) & fixPosPerTaskY{bb}>=yBins(jj) & fixPosPerTaskY{bb}<yBins(jj+1));
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
    
     subplot(3,2,3)
    histogram2Polar(relativeSaccAngleDiff1noInf,relativeSaccAmpDiff1noInf+6,.5, 'RTicks',[]); % this doesn't work because the plot won't display negative radii
    title('Model');
    
    vv = [relativeSaccAngleDiff1noInf,relativeSaccAmpDiff1noInf];
    
    ratioReturnToForewardSaccAllSaccs = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );
    disp(['Model: Ratio of Prevalence of Return To Forward Saccades = ',num2str(ratioReturnToForewardSaccAllSaccs)]);

   


%% fig. 3 equivalent


% now analyze same way as data
for bb = 1:length(d.indicesForTaskbb)
    
    % now run extract saccades and relative saccades per task
    saccsCartSample = diff(e.locSampleT{bb});
    
    [theta1,rho1] = cart2pol(saccsCartSample(:,1),saccsCartSample(:,2));
    theta1 = -1*theta1;
    
    % compute relative saccades
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
    
    % need to now wrap these around, so angles greater than 180º or less than
    % -180º end up in right place
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);
    
    % make temporal sequences plots per task
    
    % find saccades matched in amplitude less than |dr/dt|<25%
    matchedAmpIndicesF = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & relativeSaccAngleDiff1 < deg2rad(25.7143) &  relativeSaccAngleDiff1 > deg2rad(-25.7143) ); % indices here are two ahead of fixation duration vector
    matchedAmpIndicesR = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & (relativeSaccAngleDiff1 < deg2rad(-154.2857) |  relativeSaccAngleDiff1 > deg2rad(154.2857)) ); % indices here are two ahead of fixation duration vector
    
    delayVec = [0 1 2]; % corresponding to fixation 1, 2, and 3 in the sequence respectively
    
    clear fixDurMatchedAmpSampleF; clear fixDurMatchedAmpSampleR; clear meanfixDurMatchedAmpSampleF; clear stdfixDurMatchedAmpSampleF; clear nfixDurMatchedAmpSampleF; clear meanfixDurMatchedAmpSampleR; clear stdfixDurMatchedAmpSampleR; clear nfixDurMatchedAmpSampleR;
    for jj = 1:length(delayVec)
        fixDurMatchedAmpSampleF{jj} = e.durSampleT{bb}(matchedAmpIndicesF+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        fixDurMatchedAmpSampleR{jj} = e.durSampleT{bb}(matchedAmpIndicesR+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        
        % save out pertask seq data
        perTaskDurSeqModelF{bb,jj} = fixDurMatchedAmpSampleF{jj}; perTaskDurSeqModelR{bb,jj} = fixDurMatchedAmpSampleR{jj};
    end
    
    % now distances
    clear fixDurMatchedAmpSampleF; clear fixDurMatchedAmpSampleR; clear meanfixDurMatchedAmpSampleF; clear stdfixDurMatchedAmpSampleF; clear nfixDurMatchedAmpSampleF; clear meanfixDurMatchedAmpSampleR; clear stdfixDurMatchedAmpSampleR; clear nfixDurMatchedAmpSampleR;
    for jj = 1:length(delayVec)
        fixDurMatchedAmpSampleF{jj} = e.distFromFOVcenterSampleT{bb}(matchedAmpIndicesF+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        fixDurMatchedAmpSampleR{jj} = e.distFromFOVcenterSampleT{bb}(matchedAmpIndicesR+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        
        % save out pertask seq data
        perTaskDistSeqModelF{bb,jj} = fixDurMatchedAmpSampleF{jj}; perTaskDistSeqModelR{bb,jj} = fixDurMatchedAmpSampleR{jj};
    end
    
end

% now bootstrap!

% find task with maxN (IN THE DATA, not model - just to equalize total
% bootstrapped sample size in data and model for plot)
% -- copied in from Fig3.m, because we don't run this here. ***
maxNF=1149;
maxNR=4013;

% bootstrap
for bb = 1:length(d.indicesForTaskbb) % tasks
    for jj = 1:3 % # fixation
        clear bsindsF; clear bsindsF;
        
        bsindsF = randi(length(perTaskDurSeqModelF{bb,jj}),maxNF,1);
        bootstrappedDurModelF{bb,jj} = perTaskDurSeqModelF{bb,jj}(bsindsF);
        bootstrappedDistModelF{bb,jj} = perTaskDistSeqModelF{bb,jj}(bsindsF);
        
        bsindsR = randi(length(perTaskDurSeqModelR{bb,jj}),maxNR,1);
        bootstrappedDurModelR{bb,jj} = perTaskDurSeqModelR{bb,jj}(bsindsR);
        bootstrappedDistModelR{bb,jj} = perTaskDistSeqModelR{bb,jj}(bsindsR);
        
    end
end

% now pool across tasks and plot
for jj = 1:3
    bootstrappedDurModelPooledF(:,jj) = cat(2,bootstrappedDurModelF{:,jj})';
    bootstrappedDurModelPooledR(:,jj) = cat(2,bootstrappedDurModelR{:,jj})';
    bootstrappedDistModelPooledF(:,jj) = cat(2,bootstrappedDistModelF{:,jj})';
    bootstrappedDistModelPooledR(:,jj) = cat(2,bootstrappedDistModelR{:,jj})';
end


% compute SEM based on original sample size
% note that I am using an SEM error bar based on the original amount of Model, not the bootstrapped amount, which would overinflate the SEM precision
for jj = 1:3 % concat data first
    newD_durF(:,jj) = cat(2,perTaskDurSeqModelF{:,jj})';
    newD_durR(:,jj) = cat(2,perTaskDurSeqModelR{:,jj})';
    newD_distF(:,jj) = cat(2,perTaskDistSeqModelF{:,jj})';
    newD_distR(:,jj) = cat(2,perTaskDistSeqModelR{:,jj})';
end
SEM_dur_F = nanstd(newD_durF)./sqrt(size(newD_durF,1)); SEM_dur_R = nanstd(newD_durR)./sqrt(size(newD_durR,1));
SEM_dist_F = nanstd(newD_distF)./sqrt(size(newD_distF,1)); SEM_dist_R = nanstd(newD_distR)./sqrt(size(newD_distR,1));

subplot(3,2,4)
errorbar(1:3,nanmean(bootstrappedDurModelPooledF),SEM_dur_F,'ko-'); % note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
hold on; errorbar(1:3,nanmean(bootstrappedDurModelPooledR),SEM_dur_R,'ro-');
xticks([1 2 3])
xlim([.75 3.25])
ylim([.6 .75])
axis square
xlabel('Fixation # in sequence')
ylabel('Fixation duration (s)')
title('Model: Fixation Duration Sequence')

subplot(3,2,5)
errorbar(1:3,nanmean(bootstrappedDistModelPooledF),SEM_dist_F,'ko-'); % note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
hold on; errorbar(1:3,nanmean(bootstrappedDistModelPooledR),SEM_dist_R,'ro-');
xticks([1 2 3])
xlim([.75 3.25])
ylim([7 18])
axis square
xlabel('Fixation # in sequence')
ylabel('Distance from center of orbit (º)')
title('Model: Fixation Eccentricity Sequence')




%% Fig. 4  equiv

for kk = 1:length(d.fixLocationX) % loop through recordings
    fixDur2m(kk) = nanmean(fixDur2{kk}); % get mean fixation duration in each recording
end


for ii = 1:length(d.indicesForTaskbb)
    curFixDur = fixDur2m(d.indicesForTaskbb{ii});
    
    meanFixDur(ii) = nanmean(curFixDur); % median gives same result
    
    xDisp = nanstd(cat(1,fixPosX2{1, d.indicesForTaskbb{ii}})); % dispersion, std
    yDisp = nanstd(cat(1,fixPosY2{1, d.indicesForTaskbb{ii}}));
    totalDispByTask(ii) = (xDisp+yDisp)/2;
end

colorsss = hsv(9);
subplot(3,2,6)
for ii = 1:length(d.indicesForTaskbb)
    scatter(meanFixDur(ii),totalDispByTask(ii),'markerfacecolor',colorsss(ii,:)); hold on
end
[r,p] = corrcoef(meanFixDur,totalDispByTask);
text(.25,6.5,['r=' num2str(round(r(2),2)) ', p=' num2str(round(p(2),7))])
xlabel('Fixation duration (s)')
ylabel('Fixation dispersion (º)')
axis square
legend(d.taskList)


