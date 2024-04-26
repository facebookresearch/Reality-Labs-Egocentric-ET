function [] = figS8(d,e)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS8

d = d.GazeInHead;

%% DATA

saccsPooledByTask = cell(length(d.taskList),1);
for bb = 1:length(d.taskList)

    for jj = 1:length(d.indicesForTaskbb{bb})
        curInd = d.indicesForTaskbb{bb}(jj);
        saccsPooledByTask{bb} = [saccsPooledByTask{bb}; d.saccsCart{curInd}];
    end

    [theta1,rho1] = cart2pol(saccsPooledByTask{bb}(:,1),saccsPooledByTask{bb}(:,2)); % in this coordinate system, up is +90, right is 0, left, is +-180, down is -90
    theta1 = -theta1;

    % save out median sacc amp without log
    medAmpData(bb) = median(rho1);

    % Plot IOR metric
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale

    % need to now wrap these around, so angles greater than 180º or less than
    % -180º end up in right place
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);

    % estimate IOR by checking how many saccades flip sign. If below 50%
    % there may be inhibition of return
    vv = [relativeSaccAngleDiff1,relativeSaccAmpDiff1];

    ratioReturnToForewardSaccData(bb) = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );

end


%% MODEL

for bb = 1:length(d.taskList)

    saccsCartModel = diff(e.locSampleT{bb}); % compute saccades from model, and then pool 
    

    [theta1,rho1] = cart2pol(saccsCartModel(:,1),saccsCartModel(:,2)); % in this coordinate system, up is +90, right is 0, left, is +-180, down is -90
    theta1 = -theta1;

    % save out median amp without log
    medAmpModel(bb) = median(rho1);

    % Plot IOR metric
    relativeSaccAngleDiff1 = diff(theta1); % now this is in the same coordinate system as BAys and Huassain paper, where up is 0, left is -90, right +90, and down is +-180. and 0 represent a saccade in same direction
    relativeSaccAmpDiff1 = diff(log(rho1)); % difference between successive amps that have been transformed into a log scale

    % need to now wrap these around, so angles greater than 180º or less than
    % -180º end up in right place
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)>pi),1) = -pi + ( (relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)>pi,1)) - pi);
    relativeSaccAngleDiff1(find(relativeSaccAngleDiff1(:,1)<-pi),1) = pi - ( -(relativeSaccAngleDiff1(relativeSaccAngleDiff1(:,1)<-pi,1)) - pi);

    % estimate IOR by checking how many saccades flip sign. If below 50%
    % there may be inhibition of return
    vv = [relativeSaccAngleDiff1,relativeSaccAmpDiff1];

    ratioReturnToForewardSaccModel(bb) = sum( ( vv(:,1)<deg2rad(-135) | vv(:,1)>deg2rad(135) ) ) ./ sum( (vv(:,1)>deg2rad(-45) & vv(:,1)<deg2rad(45)) );

end


% plot them against each other
colorsss = hsv(9);

figure;
subplot(1,2,1);
for ii = 1:length(d.taskList)
    scatter(medAmpData(ii),medAmpModel(ii),50,'markerfacecolor',colorsss(ii,:))
    hold on;
end
xlabel('Data: saccade amplitude (º)')
ylabel('Model: saccade amplitude (º)')
xlim([4 20])
ylim([4 20])
axis square
sdsline = refline(1,0);
sdsline.Color = 'k';
legend(d.taskList)
[r,p] = corrcoef(medAmpData,medAmpModel);
text(7,5.3,['r=' num2str(r(2)) ', p=' num2str(p(2))]) 


subplot(1,2,2)
for ii = 1:length(d.taskList)
    scatter(ratioReturnToForewardSaccData(ii),ratioReturnToForewardSaccModel(ii),50,'markerfacecolor',colorsss(ii,:))
    hold on;
end
xlabel('Data: spatial IOR ratio')
ylabel('Model: spatial IOR ratio')
xlim([1.5 5])
ylim([1.5 5])
axis square
sdsline = refline(1,0);
sdsline.Color = 'k';
%legend(d.taskList)
[r,p] = corrcoef(ratioReturnToForewardSaccData,ratioReturnToForewardSaccModel);
text(2.5,2.2,['r=' num2str(r(2)) ', p=' num2str(p(2))]) 


end