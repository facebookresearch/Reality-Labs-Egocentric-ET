function [] = figS6and7(d)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS6and7

d = d.GazeInHead;

% Data: Control analysis
%% Fig. 6, panel 2
% look at whether there is temporal IOR when you match the distance from the center of the orbit for each fixation in the return and forward sequences
% if not, then it is more likely that temporal IOR is just driven by head-eye coordination + center bias + fixations longer in center of orbit (laziness)


% match the distance of the second fixation in the sequence between forward
% and return sequences and then look at their durations

clear fixDurMatchedAmp_contF; clear fixDurMatchedAmp_contR;
for kk = 1:length(d.fixLocationX)
    
    % do same as above but per recording
    [theta1,rho1] = cart2pol(d.saccsCart{kk}(:,1),d.saccsCart{kk}(:,2));
    theta1 = -1*theta1;
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
    
    % need to now wrap these around, so angles greater than 180º or less than
    % -180º end up in right place
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);
    
    % match absolute distance from center for forward and return
    
    % find saccades matched in amplitude less than |dr/dt|<25% ,
    % forward and return saccades. defined as within middle of 7 bins
    % or on end
    matchedAmpIndicesF = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & (relativeSaccAngleDiff1 < deg2rad(25.7143)) &  (relativeSaccAngleDiff1 > deg2rad(-25.7143)) ); % indices here are two ahead of fixation duration vector
    matchedAmpIndicesR = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & (relativeSaccAngleDiff1 < deg2rad(-154.2857) |  relativeSaccAngleDiff1 > deg2rad(154.2857)) ); % indices here are two ahead of fixation duration vector
    
    
    ref = median(d.distFromCenterOfFOV{kk}); % reference distance - scaled for each recording depending on distirbutions of fixation distances from center of FOV
    
    buff = 1; % this is half of interval around median; set to 1000 is equivalent to not thresholding based on distance to FOV.
    
    indsF = matchedAmpIndicesF+1;
    indsF = indsF(d.distFromCenterOfFOV{kk}((matchedAmpIndicesF+1)) < ref+buff  & d.distFromCenterOfFOV{kk}((matchedAmpIndicesF+1)) > ref-buff);
    
    indsR = matchedAmpIndicesR+1;
    indsR = indsR(d.distFromCenterOfFOV{kk}((matchedAmpIndicesR+1)) < ref+buff  & d.distFromCenterOfFOV{kk}((matchedAmpIndicesR+1)) > ref-buff);
    
    fixDurMatchedAmp_contF{kk} = d.fixDur{kk}(indsF); % but we want the duration of the second fixation in the three-fixation sequence.
    fixDurMatchedAmp_contR{kk} = d.fixDur{kk}(indsR); % but we want the duration of the second fixation in the three-fixation sequence.
end

% now bootstrap the sequence data

% pool sequences across recordings for each task
for bb = 1:length(d.indicesForTaskbb)
    for jj = 1:3
        pooledDurSeqDataF{bb} = cat(1,fixDurMatchedAmp_contF{d.indicesForTaskbb{bb}});
        pooledDurSeqDataR{bb} = cat(1,fixDurMatchedAmp_contR{d.indicesForTaskbb{bb}});
    end
end

% find task with maxN
for bb = 1:length(d.indicesForTaskbb)
    sizesssF(bb) = length(pooledDurSeqDataF{bb});
    sizesssR(bb) = length(pooledDurSeqDataR{bb});
end
[maxNF,~ ] = max(sizesssF);
[maxNR,~ ] = max(sizesssR);

for iter = 1:100
    disp(iter)
    % bootstrap
    for bb = 1:length(d.indicesForTaskbb) % tasks
        clear bsindsF; clear bsindsF;
        
        bsindsF = randi(length(pooledDurSeqDataF{bb}),maxNF,1);
        bootstrappedDurDataF{bb} = pooledDurSeqDataF{bb}(bsindsF);
        
        bsindsR = randi(length(pooledDurSeqDataR{bb}),maxNR,1);
        bootstrappedDurDataR{bb} = pooledDurSeqDataR{bb}(bsindsR);
    end
    
    % now pool across tasks and plot
    bootstrappedDurDataPooled_contF = cat(1,bootstrappedDurDataF{:});
    bootstrappedDurDataPooled_contR = cat(1,bootstrappedDurDataR{:});
    
    pPooled_cont(iter) = ptestModule(bootstrappedDurDataPooled_contF,bootstrappedDurDataPooled_contR,1000,1);
    
    meanF_cont(iter,:) = nanmean(bootstrappedDurDataPooled_contF);
    semF_cont(iter,:) = nanstd(bootstrappedDurDataPooled_contF)./sqrt(length(bootstrappedDurDataPooled_contF));
    
    meanR_cont(iter,:) = nanmean(bootstrappedDurDataPooled_contR);
    semR_cont(iter,:) = nanstd(bootstrappedDurDataPooled_contR)./sqrt(length(bootstrappedDurDataPooled_contR));
    
end

saveOutSizeF = length(bootstrappedDurDataPooled_contF); % for matching raw to control analysis sample size
saveOutSizeR = length(bootstrappedDurDataPooled_contR);

figure;
subplot(1,2,2)
errorbar(1, mean(meanF_cont),mean(semF_cont),'ko-' )
hold on; errorbar(1.4, mean(meanR_cont), mean(semR_cont),'ro-' )
xlim([.75 1.7])
ylim([.41 .55])
xticks([1,1.4])
xticklabels({'forward','return'})
legend('Forward','Return')
xlabel('Fixation # in sequence')
ylabel('Fixation duration (s)')
title('Distance Matched')



%% Fig. 6- Panel 1 : normal data but sample size matched to original data. 

clear pooledDurSeqDataF; clear pooledDurSeqDataR;
for kk = 1:length(d.fixLocationX)
    
    % do same as above but per recording
    [theta1,rho1] = cart2pol(d.saccsCart{kk}(:,1),d.saccsCart{kk}(:,2));
    theta1 = -1*theta1;
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
    
    % need to now wrap these around, so angles greater than 180º or less than
    % -180º end up in right place
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);
    
    % match absolute distance from center for forward and return
    
    % find saccades matched in amplitude less than |dr/dt|<25% ,
    % forward and return saccades. defined as within middle of 7 bins
    % or on end
    matchedAmpIndicesF = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & (relativeSaccAngleDiff1 < deg2rad(25.7143)) &  (relativeSaccAngleDiff1 > deg2rad(-25.7143)) ); % indices here are two ahead of fixation duration vector
    matchedAmpIndicesR = find(abs(relativeSaccAmpDiff1)<=quantile(abs(relativeSaccAmpDiff1),.25) & (relativeSaccAngleDiff1 < deg2rad(-154.2857) |  relativeSaccAngleDiff1 > deg2rad(154.2857)) ); % indices here are two ahead of fixation duration vector
    
    
    ref = median(d.distFromCenterOfFOV{kk}); % reference distance - scaled for each recording depending on distirbutions of fixation distances from center of FOV
    
    buff = 1000; % this is half of interval around median; set to 1000 is equivalent to not thresholding based on distance to FOV.
    
    indsF = matchedAmpIndicesF+1;
    indsF = indsF(d.distFromCenterOfFOV{kk}((matchedAmpIndicesF+1)) < ref+buff  & d.distFromCenterOfFOV{kk}((matchedAmpIndicesF+1)) > ref-buff);
    
    indsR = matchedAmpIndicesR+1;
    indsR = indsR(d.distFromCenterOfFOV{kk}((matchedAmpIndicesR+1)) < ref+buff  & d.distFromCenterOfFOV{kk}((matchedAmpIndicesR+1)) > ref-buff);
    
    fixDurMatchedAmpF{kk} = d.fixDur{kk}(indsF); % but we want the duration of the second fixation in the three-fixation sequence.
    fixDurMatchedAmpR{kk} = d.fixDur{kk}(indsR); % but we want the duration of the second fixation in the three-fixation sequence.
end

% now bootstrap the sequence data

% pool sequences across recordings for each task
for bb = 1:length(d.indicesForTaskbb)
    for jj = 1:3
        pooledDurSeqDataF{bb} = cat(1,fixDurMatchedAmpF{d.indicesForTaskbb{bb}});
        pooledDurSeqDataR{bb} = cat(1,fixDurMatchedAmpR{d.indicesForTaskbb{bb}});
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
    clear bsindsF; clear bsindsF;
    
    bsindsF = randi(length(pooledDurSeqDataF{bb}),maxNF,1);
    bootstrappedDurDataF{bb} = pooledDurSeqDataF{bb}(bsindsF);
    
    bsindsR = randi(length(pooledDurSeqDataR{bb}),maxNR,1);
    bootstrappedDurDataR{bb} = pooledDurSeqDataR{bb}(bsindsR);
end

% now pool across tasks and plot
bootstrappedDurDataPooledF = cat(1,bootstrappedDurDataF{:});
bootstrappedDurDataPooledR = cat(1,bootstrappedDurDataR{:});

% pooled fig
% uncomment - for when you add 1000 to the bounds above to include
% datapoints at all distances from center of FOV. this loop is to match
% data sample size for distance restricted and non restricted case. and see
% how sampling error effects p-value meta distribution.
for pIter = 1:100
    disp(pIter)
    clear reducedF; clear reducedR;
    reducedF = bootstrappedDurDataPooledF(randi(length(bootstrappedDurDataPooledF),saveOutSizeF,1));
    reducedR = bootstrappedDurDataPooledR(randi(length(bootstrappedDurDataPooledR),saveOutSizeR,1));

    meanF = mean(reducedF); semF(pIter) = std(reducedF)./sqrt(length(reducedF));
    meanR = mean(reducedR); semR(pIter) = std(reducedR)./sqrt(length(reducedR));
    
    pPooled(pIter) = ptestModule(reducedF,reducedR,1000,1);
end

subplot(1,2,1)
errorbar(1, mean(meanF),mean(semF),'ko-' )
hold on; errorbar(1.4, mean(meanR), mean(semR),'ro-' )
xlim([.75 1.7])
ylim([.41 .55])
xticks([1,1.4])
xticklabels({'forward','return'})
legend('Forward','Return')
xlabel('Fixation # in sequence')
ylabel('Fixation duration (s)')
title('Raw (sample size matched)')

%{
subplot(1,2,1)
errorbar(1, mean(bootstrappedDurDataPooledF),std(bootstrappedDurDataPooledF)./sqrt(length(bootstrappedDurDataPooledF)),'ko-' )
hold on; errorbar(1.4, mean(bootstrappedDurDataPooledR), std(bootstrappedDurDataPooledR)./sqrt(length(bootstrappedDurDataPooledR)),'ro-' )
xlim([.75 1.7])
xticks([1,1.4])
xticklabels({'forward','return'})
legend('Forward','Return')
xlabel('Fixation # in sequence')
ylabel('Fixation duration (s)')
%}




%% Fig. 7 - per-task figs (based on control analysis in Fig.S6a)

clear p;
% per task figs
for ii = 1:length(d.indicesForTaskbb)
    
    perTaskControlDursF{ii} = cat(1,fixDurMatchedAmp_contF{d.indicesForTaskbb{ii}});
    perTaskControlDursR{ii} = cat(1,fixDurMatchedAmp_contR{d.indicesForTaskbb{ii}});
    
    meanPerTaskControlDursF(ii,:) = nanmean(perTaskControlDursF{ii});
    stdPerTaskControlDursF(ii,:) = nanstd(perTaskControlDursF{ii});
    nPerTaskControlDursF(ii,:) = length(perTaskControlDursF{ii});
    
    meanPerTaskControlDursR(ii,:) = nanmean(perTaskControlDursR{ii});
    stdPerTaskControlDursR(ii,:) = nanstd(perTaskControlDursR{ii});
    nPerTaskControlDursR(ii,:) =  length(perTaskControlDursR{ii});
    
    
    % do permutation tests
    p(ii) = ptestModule( perTaskControlDursR{ii},perTaskControlDursF{ii},1000,0);
end


% per task IOR sequence
figure;
for ii = 1:length(d.indicesForTaskbb)
    
    subplot(3,3,ii); errorbar(1, meanPerTaskControlDursF(ii,:),stdPerTaskControlDursF(ii,:)./sqrt(nPerTaskControlDursF(ii,:)),'ko-' )
    hold on; errorbar(1.4, meanPerTaskControlDursR(ii,:),stdPerTaskControlDursR(ii,:)./sqrt(nPerTaskControlDursR(ii,:)),'ro-' )
    xticks([1,1.4])
    xticklabels({'forward','return'})
    xlim([.75 1.7])
    axis square
    title(d.taskList{ii})
    
    if ii == 7
        legend('Forward','Return')
        %ylim([.4 .6])
        xlabel('Fixation # in sequence')
        ylabel('Fixation duration (s)')
    end
    
end

suptitle('Distance Matched')



end