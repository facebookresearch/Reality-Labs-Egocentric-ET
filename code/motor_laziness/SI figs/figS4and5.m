function [] = figS4and5(d,e)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS4and5

d = d.GazeInHead;


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
        
        % now distances
        clear fixDurMatchedAmpSampleF; clear fixDurMatchedAmpSampleR; clear meanfixDurMatchedAmpSampleF; clear stdfixDurMatchedAmpSampleF; clear nfixDurMatchedAmpSampleF; clear meanfixDurMatchedAmpSampleR; clear stdfixDurMatchedAmpSampleR; clear nfixDurMatchedAmpSampleR;
        for jj = 1:length(delayVec)
            fixDurMatchedAmpSampleF{jj} = d.distFromCenterOfFOV{kk}(matchedAmpIndicesF+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
            fixDurMatchedAmpSampleR{jj} = d.distFromCenterOfFOV{kk}(matchedAmpIndicesR+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
            
            % save out per-recording seq data
            perRecordingDistSeqDataF{kk,jj} = fixDurMatchedAmpSampleF{jj}; perRecordingDistSeqDataR{kk,jj} = fixDurMatchedAmpSampleR{jj};
        end
        
    end
    
end

% pool sequences across recordings for each task
for bb = 1:length(d.indicesForTaskbb)
    for jj = 1:3
        pooledDurSeqDataF{bb}(:,jj) = cat(1,perRecordingDurSeqDataF{d.indicesForTaskbb{bb},jj});
        pooledDurSeqDataR{bb}(:,jj) = cat(1,perRecordingDurSeqDataR{d.indicesForTaskbb{bb},jj});
        pooledDistSeqDataF{bb}(:,jj) = cat(1,perRecordingDistSeqDataF{d.indicesForTaskbb{bb},jj});
        pooledDistSeqDataR{bb}(:,jj) = cat(1,perRecordingDistSeqDataR{d.indicesForTaskbb{bb},jj});
    end
end

% per task IOR sequence
figure;
for ii = 1:length(d.indicesForTaskbb)
    
    subplot(3,3,ii); errorbar(1:3,nanmean(pooledDurSeqDataF{ii}),nanstd(pooledDurSeqDataF{ii})./sqrt(size(pooledDurSeqDataF{ii},1)),'ko-');
    hold on;  errorbar(1:3,nanmean(pooledDurSeqDataR{ii}),nanstd(pooledDurSeqDataR{ii})./sqrt(size(pooledDurSeqDataR{ii},1)),'ro-');
    xticks([1 2 3])
    xlim([.75 3.25])
    axis square
    title(d.taskList{ii})
    
    if ii == 7
        %legend('Forward','Return')
        %ylim([.4 .6])
        xlabel('Fixation # in sequence')
        ylabel('Fixation duration (s)')
    end
end
suptitle('Data')

figure;
for ii = 1:length(d.indicesForTaskbb)
    
    subplot(3,3,ii); errorbar(1:3,nanmean(pooledDistSeqDataF{ii}),nanstd(pooledDistSeqDataF{ii})./sqrt(size(pooledDistSeqDataF{ii},1)),'ko-');
    hold on;  errorbar(1:3,nanmean(pooledDistSeqDataR{ii}),nanstd(pooledDistSeqDataR{ii})./sqrt(size(pooledDistSeqDataR{ii},1)),'ro-');
    xticks([1 2 3])
    xlim([.75 3.25])
    axis square
    title(d.taskList{ii})
    
    if ii == 7
        %legend('Forward','Return')
        %ylim([.4 .6])
        xlabel('Fixation # in sequence')
        ylabel('Distance from center of orbit (º)')
    end
end

suptitle('Data')



%% MODEL: temporal IOR sequence summary fig

% load in simulated fixations and process them in the same way


% first remove excess simulated data per task. 
for bb = 1:length(d.indicesForTaskbb)
    perTaskFixDur = cat(1,d.fixDur{d.indicesForTaskbb{bb}});
    
    howMuchDataToRemove = length(e.durSampleT{bb})-length(perTaskFixDur);
    whichSamplesToRemove = randi(length(e.durSampleT{bb}),howMuchDataToRemove,1);
    
    % remove model data randomly to match original data sample size
    e.locSampleT2{bb}=e.locSampleT{bb};
    e.locSampleT2{bb}(whichSamplesToRemove,:) = [];
    
    e.durSampleT2{bb}=e.durSampleT{bb};
    e.durSampleT2{bb}(whichSamplesToRemove) = [];
    
    e.distFromFOVcenterSampleT2{bb}=e.distFromFOVcenterSampleT{bb};
    e.distFromFOVcenterSampleT2{bb}(whichSamplesToRemove) = [];
end


for bb = 1:length(d.indicesForTaskbb)
    
    % now run extract saccades and relative saccades per task
    saccsCartSample = diff(e.locSampleT2{bb});
    
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
        fixDurMatchedAmpSampleF{jj} = e.durSampleT2{bb}(matchedAmpIndicesF+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        fixDurMatchedAmpSampleR{jj} = e.durSampleT2{bb}(matchedAmpIndicesR+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        
        meanfixDurMatchedAmpSampleF(jj) = nanmean(fixDurMatchedAmpSampleF{jj});
        stdfixDurMatchedAmpSampleF(jj) = nanstd(fixDurMatchedAmpSampleF{jj});
        nfixDurMatchedAmpSampleF(jj) = length(fixDurMatchedAmpSampleF{jj});
        
        meanfixDurMatchedAmpSampleR(jj) = nanmean(fixDurMatchedAmpSampleR{jj});
        stdfixDurMatchedAmpSampleR(jj) = nanstd(fixDurMatchedAmpSampleR{jj});
        nfixDurMatchedAmpSampleR(jj) = length(fixDurMatchedAmpSampleR{jj});
    end
    
    % save out pertask seq data
    meanfixDurMatchedAmpSamplePerTaskF(bb,:) = meanfixDurMatchedAmpSampleF; meanfixDurMatchedAmpSamplePerTaskR(bb,:) = meanfixDurMatchedAmpSampleR;
    stdfixDurMatchedAmpSamplePerTaskF(bb,:) = stdfixDurMatchedAmpSampleF; stdfixDurMatchedAmpSamplePerTaskR(bb,:) = stdfixDurMatchedAmpSampleR;
    nfixDurMatchedAmpSamplePerTaskF(bb,:) = nfixDurMatchedAmpSampleF; nfixDurMatchedAmpSamplePerTaskR(bb,:) = nfixDurMatchedAmpSampleR;
    
    clear fixDurMatchedAmpSampleF; clear fixDurMatchedAmpSampleR; clear meanfixDurMatchedAmpSampleF; clear stdfixDurMatchedAmpSampleF; clear nfixDurMatchedAmpSampleF; clear meanfixDurMatchedAmpSampleR; clear stdfixDurMatchedAmpSampleR; clear nfixDurMatchedAmpSampleR;
    for jj = 1:length(delayVec)
        fixDurMatchedAmpSampleF{jj} = e.distFromFOVcenterSampleT2{bb}(matchedAmpIndicesF+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        fixDurMatchedAmpSampleR{jj} = e.distFromFOVcenterSampleT2{bb}(matchedAmpIndicesR+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        
        meanfixDurMatchedAmpSampleF(jj) = nanmean(fixDurMatchedAmpSampleF{jj});
        stdfixDurMatchedAmpSampleF(jj) = nanstd(fixDurMatchedAmpSampleF{jj});
        nfixDurMatchedAmpSampleF(jj) = length(fixDurMatchedAmpSampleF{jj});
        
        meanfixDurMatchedAmpSampleR(jj) = nanmean(fixDurMatchedAmpSampleR{jj});
        stdfixDurMatchedAmpSampleR(jj) = nanstd(fixDurMatchedAmpSampleR{jj});
        nfixDurMatchedAmpSampleR(jj) = length(fixDurMatchedAmpSampleR{jj});
    end
    
    %save out pertask seq data distFOV
    meandistFOVMatchedAmpSamplePerTaskF(bb,:) = meanfixDurMatchedAmpSampleF; meandistFOVMatchedAmpSamplePerTaskR(bb,:) = meanfixDurMatchedAmpSampleR;
    stddistFOVMatchedAmpSamplePerTaskF(bb,:) = stdfixDurMatchedAmpSampleF; stddistFOVMatchedAmpSamplePerTaskR(bb,:) = stdfixDurMatchedAmpSampleR;
    ndistFOVMatchedAmpSamplePerTaskF(bb,:) = nfixDurMatchedAmpSampleF; ndistFOVMatchedAmpSamplePerTaskR(bb,:) = nfixDurMatchedAmpSampleR;
    
end


% plot per task IOR sequences
figure;
for ii = 1:length(d.indicesForTaskbb)
    
    subplot(3,3,ii); errorbar(1:3,meanfixDurMatchedAmpSamplePerTaskF(ii,:),stdfixDurMatchedAmpSamplePerTaskF(ii,:)./sqrt(nfixDurMatchedAmpSamplePerTaskF(ii,:)),'ko-');
    hold on; errorbar(1:3,meanfixDurMatchedAmpSamplePerTaskR(ii,:),stdfixDurMatchedAmpSamplePerTaskR(ii,:)./sqrt(nfixDurMatchedAmpSamplePerTaskR(ii,:)),'ro-');
    xticks([1 2 3])
    xlim([.75 3.25])
    axis square
    title(d.taskList{ii})
    
    if ii == 7
        %legend('Foreward','Return')
        %ylim([.4 .6])
        xlabel('Fixation # in sequence')
        ylabel('Fixation duration (s)')
    end
    
end
suptitle('Model')

figure;
for ii = 1:length(d.indicesForTaskbb)
    
    subplot(3,3,ii); errorbar(1:3,meandistFOVMatchedAmpSamplePerTaskF(ii,:),stddistFOVMatchedAmpSamplePerTaskF(ii,:)./sqrt(ndistFOVMatchedAmpSamplePerTaskF(ii,:)),'ko-');
    hold on; errorbar(1:3,meandistFOVMatchedAmpSamplePerTaskR(ii,:),stddistFOVMatchedAmpSamplePerTaskR(ii,:)./sqrt(ndistFOVMatchedAmpSamplePerTaskR(ii,:)),'ro-');
    xticks([1 2 3])
    xlim([.75 3.25])
    axis square
    title(d.taskList{ii})
    
    if ii == 7
        %legend('Foreward','Return')
        %ylim([.4 .6])
        xlabel('Fixation # in sequence')
        ylabel('Distance from center of orbit (º)')
    end
    
end
suptitle('Model')


end
