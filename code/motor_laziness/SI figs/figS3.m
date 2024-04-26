function [] = figS3(d,e)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS3

d = d.GazeInHead;

%% DATA: plot per-task joint distributions

figure;
saccsPooledByTask = cell(length(d.taskList),1);
for bb = 1:length(d.taskList)
    
    for jj = 1:length(d.indicesForTaskbb{bb})
        curInd = d.indicesForTaskbb{bb}(jj);
        saccsPooledByTask{bb} = [saccsPooledByTask{bb}; d.saccsCart{curInd}];
    end
    
    [theta1,rho1] = cart2pol(saccsPooledByTask{bb}(:,1),saccsPooledByTask{bb}(:,2)); % in this coordinate system, up is +90, right is 0, left, is +-180, down is -90
    theta1 = -theta1;
    
    % Plot IOR metric
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale
    
    % need to now wrap these around, so angles greater than 180º or less than
    % -180º end up in right place
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);
    
     
    vv = [relativeSaccAngleDiff1,relativeSaccAmpDiff1];
    percentFlipSign(bb) = sum(vv(:,1)>pi/2 | vv(:,1)<-pi/2)./length(vv(:,1));
    ratioReturnToForewardSaccD(bb) = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );
    
    % plot
    subplot(3,3,bb);
    histogram2Polar(relativeSaccAngleDiff1,relativeSaccAmpDiff1+6,.5, 'RTicks',[]); % this doesn't work because the plot won't display negative radii
    colorbar off
    title(d.taskList{bb})
   
    
end

suptitle('Data')



%% MODEL: same thing

figure;
for bb = 1:length(d.taskList)
    
    % compute saccades from synthetic data
    % now run analyses per task
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
    
    % validation check: IOR ratio
    clear vv;
    vv = [relativeSaccAngleDiff1,relativeSaccAmpDiff1];
    ratioReturnToForewardSaccM(bb) = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );
    
    
    % must remove inf and inf values from log transform ...
    relativeSaccAmpDiff1noInf = relativeSaccAmpDiff1;
    relativeSaccAngleDiff1noInf = relativeSaccAngleDiff1;
    infInds = (relativeSaccAmpDiff1noInf==inf | relativeSaccAmpDiff1noInf==-inf);
    relativeSaccAmpDiff1noInf(infInds) = [];
    relativeSaccAngleDiff1noInf(infInds) = [];
    
    % draw relative saccade heatmaps
    subplot(3,3,bb)
    histogram2Polar(relativeSaccAngleDiff1noInf,relativeSaccAmpDiff1noInf+6,.5, 'RTicks',[]); % this doesn't work because the plot won't display negative radii
    colorbar off
    title(d.taskList{bb})
end

suptitle('Model')



end