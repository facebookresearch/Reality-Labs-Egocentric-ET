function [] = figS10(d)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS10

d = d.GazeInSpace; % Now using Gaze-in-space data

%% DATA

clear saccsPooledByTaskBootstrap

% find the N from the task with the highest N, that will be the number of
% fixations we draw from each task to equalize data across tasks
for ii = 1:length(d.taskList); xvx(ii) = size(d.saccsPooledByTask{ii},1); end
maxN = max(xvx);

for ii = 1:length(d.taskList)
    [theta1,rho1] = cart2pol(d.saccsPooledByTask{ii}(:,1),d.saccsPooledByTask{ii}(:,2)); % in this coordinate system, up is +90, right is 0, left, is +-180, down is -90
    theta1 = -theta1;
    
    medAmpData(ii) = nanmedian(rho1);
    semAmpData(ii) = nanstd(rho1)./sqrt(length(rho1));

    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
    
      % need to now wrap these around, so angles greater than 180º or less than
    % -180º end up in right place
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);

    % estimate IOR by checking how many saccades flip sign. If below 50%
    % there may be inhibition of return
    vv = [relativeSaccAngleDiff1,relativeSaccAmpDiff1];

    ratioReturnToForewardSaccData(ii) = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );

    both = [relativeSaccAngleDiff1 relativeSaccAmpDiff1];
    
    bs_inds = randi(size(both,1),maxN,1);
    saccsPooledByTaskBootstrap{ii} = both(bs_inds,:);
end


allSaccDiffs = cat(1,saccsPooledByTaskBootstrap{:});

relativeSaccAngleDiff1 = allSaccDiffs(:,1);
relativeSaccAmpDiff1 = allSaccDiffs(:,2);


% plot sacc amplitudes across tasks
subplot(3,2,1); bar(medAmpData,'k')
hold on;
errorbar(1:9,medAmpData,semAmpData,'k.')
xticklabels(d.taskList)
xtickangle(45) 
ylabel('Median saccade amplitude (º)')
xlabel('Task')

% plot spatial IOR ratio across tasks
subplot(3,2,2); bar(ratioReturnToForewardSaccData,'k')
xticklabels(d.taskList)
xtickangle(45) 
ylabel('Spatial IOR ratio')
xlabel('Task')

%% Plot pooled tasks relative saccade joint distribution

% need to now wrap these around, so angles greater than 180º or less than
% -180º end up in right place
relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);

 subplot(3,2,3);
histogram2Polar(relativeSaccAngleDiff1,relativeSaccAmpDiff1+6,.5); % this doesn't work because the plot won't display negative radii
vv = [relativeSaccAngleDiff1,relativeSaccAmpDiff1];
title('Data')
colorbar

ratioReturnToForewardSaccAllSaccs = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );

disp(['Data: Ratio of Prevalence of Return To Forward Saccades = ',num2str(ratioReturnToForewardSaccAllSaccs)]);



%% DATA: temporal IOR sequence summary fig

delayVec = [0 1 2]; % corresponding to fixation 1, 2, and 3 in the sequence respectively

for jj = 1:length(delayVec)
    
    for kk = 1:length(d.fixDur)
        
        % do same as above but per recording
        [theta1,rho1] = cart2pol(d.saccsCart{kk}(:,1),d.saccsCart{kk}(:,2));
        theta1 = -1*theta1;
        relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
        relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
        
        % need to now wrap these around, so angles greater than 180º or less than
        % -180º end up in right place
        relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
        relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);
             
        % find saccades matched in amplitude less than |dr/dt|<25%
        matchedAmpIndicesF = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & relativeSaccAngleDiff1 < deg2rad(25.7143) &  relativeSaccAngleDiff1 > deg2rad(-25.7143) ); % indices here are two ahead of fixation duration vector
        matchedAmpIndicesR = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & (relativeSaccAngleDiff1 < deg2rad(-154.2857) |  relativeSaccAngleDiff1 > deg2rad(154.2857)) ); % indices here are two ahead of fixation duration vector
        
        delayVec = [0 1 2]; % corresponding to fixation 1, 2, and 3 in the sequence respectively
        
        clear fixDurMatchedAmpSampleF; clear fixDurMatchedAmpSampleR; clear meanfixDurMatchedAmpSampleF; clear stdfixDurMatchedAmpSampleF; clear nfixDurMatchedAmpSampleF; clear meanfixDurMatchedAmpSampleR; clear stdfixDurMatchedAmpSampleR; clear nfixDurMatchedAmpSampleR;
        for jj = 1:length(delayVec)
            fixDurMatchedAmpSampleF{jj} = d.fixDur{kk}(matchedAmpIndicesF+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
            fixDurMatchedAmpSampleR{jj} = d.fixDur{kk}(matchedAmpIndicesR+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
            
            % save out per-recording seq data
            perRecordingDurSeqDataF{kk,jj} = fixDurMatchedAmpSampleF{jj}; perRecordingDurSeqDataR{kk,jj} = fixDurMatchedAmpSampleR{jj};
        end
        
           
    end
    
end

% pool sequences across recordings for each task
for bb = 1:length(d.indicesForTaskbb)
    for jj = 1:3
        pooledDurSeqDataF{bb}(:,jj) = cat(1,perRecordingDurSeqDataF{d.indicesForTaskbb{bb},jj});
        pooledDurSeqDataR{bb}(:,jj) = cat(1,perRecordingDurSeqDataR{d.indicesForTaskbb{bb},jj});
    end
end

% find task with maxN
for bb = 1:length(d.indicesForTaskbb)
    sizesssF(bb) = length(pooledDurSeqDataF{bb});
    sizesssR(bb) = length(pooledDurSeqDataR{bb});
end
[maxNF,~ ] = max(sizesssF);
[maxNR,~ ] = max(sizesssR);

% bootstrap
for bb = 1:length(d.indicesForTaskbb) % tasks
    for jj = 1:3 % # fixation
        clear bsindsF; clear bsindsF;
        
        bsindsF = randi(length(pooledDurSeqDataF{bb}),maxNF,1);
        bootstrappedDurDataF{bb,jj} = pooledDurSeqDataF{bb}(bsindsF,jj);
        
        bsindsR = randi(length(pooledDurSeqDataR{bb}),maxNR,1);
        bootstrappedDurDataR{bb,jj} = pooledDurSeqDataR{bb}(bsindsR,jj);
        
    end
end

% now pool across tasks and plot
for jj = 1:3
    bootstrappedDurDataPooledF(:,jj) = cat(1,bootstrappedDurDataF{:,jj});
    bootstrappedDurDataPooledR(:,jj) = cat(1,bootstrappedDurDataR{:,jj});
end

% compute SEM based on original sample size
% note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
SEM_dur_F = std(cat(1,pooledDurSeqDataF{:}))./sqrt(size(cat(1,pooledDurSeqDataF{:}),1)); SEM_dur_R = std(cat(1,pooledDurSeqDataR{:}))./sqrt(size(cat(1,pooledDurSeqDataR{:}),1));


subplot(3,2,4)
errorbar(1:3,nanmean(bootstrappedDurDataPooledF),SEM_dur_F,'ko-'); % note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
hold on; errorbar(1:3,nanmean(bootstrappedDurDataPooledR),SEM_dur_R,'ro-');
xticks([1 2 3])
xlim([.75 3.25])
ylim([.38 .65])
axis square
xlabel('Fixation # in sequence')
ylabel('Fixation duration (s)')


% Pvals before bootstrapping
for bb = 1:9
    pPerTask(bb)= ptestModule(pooledDurSeqDataR{bb}(:,2),pooledDurSeqDataF{bb}(:,2),1000,1);
end


%% gaze-in-space fixation location and duration maps

for bb = 1:length(d.indicesForTaskbb)
    
    fixPosPerTaskX{bb} = cat(1,d.fixLocationX{d.indicesForTaskbb{bb}});
    fixPosPerTaskY{bb} = cat(1,d.fixLocationY{d.indicesForTaskbb{bb}});
    fixDurPerTask{bb} = cat(1,d.fixDur{d.indicesForTaskbb{bb}});
end

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

% pool data across tasks
posAlllll = cat(1,fixPosPerTaskBootstrap{:});
fixAllll= cat(1,fixDurPerTaskBootstrap{:});

% spatial binning
xBins = linspace(-80,80,60);
yBins = linspace(-80,80,60);

% generate fixation position and duration maps
for ii = 1:length(xBins)-1
    for jj = 1:length(yBins)-1
        justDataInBin = fixAllll(posAlllll(:,1)>=xBins(ii) & posAlllll(:,1)<xBins(ii+1) & posAlllll(:,2)>=yBins(jj) & posAlllll(:,2)<yBins(jj+1));
        meanmap1(jj,ii) = nanmean(justDataInBin);
        semMap(jj,ii) = nanstd(justDataInBin)./sqrt(length(justDataInBin));
        nPointsmap1(jj,ii) = length(justDataInBin);
    end
end


subplot(3,2,5);
imagesc(xBins,yBins,nPointsmap1)
axis equal; hline(0,'-k'); vline(0,'-k');
set(gca,'YDir','normal')
colormap(flipud(parula))
xlim([-80 80])
ylim([-80 80])
colorbar

subplot(3,2,6);
imagesc(xBins,yBins,meanmap1,'AlphaData',nPointsmap1>=22 ) % 22 obeys the 2.2 data:bootstrappedData ratio, equivalent to removing 10 datapoints
axis equal; hline(0,'-k'); vline(0,'-k');
set(gca,'YDir','normal')
colormap(flipud(parula))
caxis([.1 .8])
xlim([-80 80])
ylim([-80 80])
colorbar



%%


end