function [] = fig3(d,e)
% Copyright (c) Meta Platforms, Inc. and affiliates.
%fig3

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
        bootstrappedDistDataF{bb,jj} = pooledDistSeqDataF{bb}(bsindsF,jj);
        
        bsindsR = randi(length(pooledDurSeqDataR{bb}),maxNR,1);
        bootstrappedDurDataR{bb,jj} = pooledDurSeqDataR{bb}(bsindsR,jj);
        bootstrappedDistDataR{bb,jj} = pooledDistSeqDataR{bb}(bsindsR,jj);
        
    end
end

% now pool across tasks and plot
for jj = 1:3
    bootstrappedDurDataPooledF(:,jj) = cat(1,bootstrappedDurDataF{:,jj});
    bootstrappedDurDataPooledR(:,jj) = cat(1,bootstrappedDurDataR{:,jj});
    bootstrappedDistDataPooledF(:,jj) = cat(1,bootstrappedDistDataF{:,jj});
    bootstrappedDistDataPooledR(:,jj) = cat(1,bootstrappedDistDataR{:,jj});
end

% compute SEM based on original sample size
% note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
SEM_dur_F = std(cat(1,pooledDurSeqDataF{:}))./sqrt(size(cat(1,pooledDurSeqDataF{:}),1)); SEM_dur_R = std(cat(1,pooledDurSeqDataR{:}))./sqrt(size(cat(1,pooledDurSeqDataR{:}),1));
SEM_dist_F = std(cat(1,pooledDistSeqDataF{:}))./sqrt(size(cat(1,pooledDistSeqDataF{:}),1)); SEM_dist_R = std(cat(1,pooledDistSeqDataR{:}))./sqrt(size(cat(1,pooledDistSeqDataR{:}),1));


figure(3);
subplot(2,2,1)
errorbar(1:3,nanmean(bootstrappedDurDataPooledF),SEM_dur_F,'ko-'); % note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
hold on; errorbar(1:3,nanmean(bootstrappedDurDataPooledR),SEM_dur_R,'ro-');
xticks([1 2 3])
xlim([.75 3.25])
ylim([.38 .65])
axis square
xlabel('Fixation # in sequence')
ylabel('Fixation duration (s)')
title('Data: Fixation Duration Sequence')

figure(3);
subplot(2,2,3)
errorbar(1:3,nanmean(bootstrappedDistDataPooledF),SEM_dist_F,'ko-'); % note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
hold on; errorbar(1:3,nanmean(bootstrappedDistDataPooledR),SEM_dist_R,'ro-');
xticks([1 2 3])
xlim([.75 3.25])
ylim([7 18])
axis square
xlabel('Fixation # in sequence')
ylabel('Distance from center of orbit (º)')
title('Data: Fixation Eccentricity Sequence')



%% MODEL: input saved out data from model, but first remove excess data per task, b/c were going to bootstrap output data same way we do with data

clear perTaskDistSeqModelF; clear perTaskDistSeqModelR; clear perTaskDurSeqModelF; clear perTaskDurSeqModelR;
clear pooledDurSeqModelF; clear pooledDurSeqModelR; clear pooledDistSeqModelF; clear pooledDistSeqModelR;

% remove excess simulated data per task. 
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

% now analyze same way as data
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
        
        % save out pertask seq data
        perTaskDurSeqModelF{bb,jj} = fixDurMatchedAmpSampleF{jj}; perTaskDurSeqModelR{bb,jj} = fixDurMatchedAmpSampleR{jj};
    end
    
    % now distances
    clear fixDurMatchedAmpSampleF; clear fixDurMatchedAmpSampleR; clear meanfixDurMatchedAmpSampleF; clear stdfixDurMatchedAmpSampleF; clear nfixDurMatchedAmpSampleF; clear meanfixDurMatchedAmpSampleR; clear stdfixDurMatchedAmpSampleR; clear nfixDurMatchedAmpSampleR;
    for jj = 1:length(delayVec)
        fixDurMatchedAmpSampleF{jj} = e.distFromFOVcenterSampleT2{bb}(matchedAmpIndicesF+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        fixDurMatchedAmpSampleR{jj} = e.distFromFOVcenterSampleT2{bb}(matchedAmpIndicesR+delayVec(jj)); % but we want the duration of the second fixation in the three-fixation sequence.
        
        % save out pertask seq data
        perTaskDistSeqModelF{bb,jj} = fixDurMatchedAmpSampleF{jj}; perTaskDistSeqModelR{bb,jj} = fixDurMatchedAmpSampleR{jj};
    end
    
end

% now bootstrap!

% find task with maxN (IN THE DATA, not model - just to equalize total
% bootstrapped sample size in data and model for plot)
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

figure(3);
subplot(2,2,2)
errorbar(1:3,nanmean(bootstrappedDurModelPooledF),SEM_dur_F,'ko-'); % note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
hold on; errorbar(1:3,nanmean(bootstrappedDurModelPooledR),SEM_dur_R,'ro-');
xticks([1 2 3])
xlim([.75 3.25])
ylim([.38 .65])
axis square
xlabel('Fixation # in sequence')
ylabel('Fixation duration (s)')
title('Model: Fixation Duration Sequence')

figure(3);
subplot(2,2,4)
errorbar(1:3,nanmean(bootstrappedDistModelPooledF),SEM_dist_F,'ko-'); % note that I am using an SEM error bar based on the original amount of data, not the bootstrapped amount, which would overinflate the SEM precision
hold on; errorbar(1:3,nanmean(bootstrappedDistModelPooledR),SEM_dist_R,'ro-');
xticks([1 2 3])
xlim([.75 3.25])
ylim([7 18])
axis square
xlabel('Fixation # in sequence')
ylabel('Distance from center of orbit (º)')
title('Model: Fixation Eccentricity Sequence')




end
