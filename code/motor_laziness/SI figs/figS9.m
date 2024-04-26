function [] = figS9(d)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS9

d = d.GazeInHead;

%% DATA

figure;
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
    semAmpData(bb) = std(rho1)./sqrt(length(rho1));
    

    subplot(3,3,bb)
    H1 = polarhistogram(theta1,40,'Normalization','pdf');
    title(d.taskList{bb})

end

figure; bar(medAmpData,'k')
hold on;
errorbar(1:9,medAmpData,semAmpData,'k.')
xticklabels(d.taskList)
xtickangle(45) 
ylabel('Median saccade amplitude (º)')
xlabel('Task')

end

