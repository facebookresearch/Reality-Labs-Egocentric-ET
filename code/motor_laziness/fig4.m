function [] = fig4(d)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% fig 4

d = d.GazeInHead;

%% DATA: Create correlation plot between avg. fixation duration and dispersion across tasks

for kk = 1:length(d.fixLocationX) % loop through recordings
    fixDur2(kk) = nanmean(d.fixDur{kk}); % get mean fixation duration in each recording
end


for ii = 1:length(d.indicesForTaskbb)
    curFixDur = fixDur2(d.indicesForTaskbb{ii});
    
    meanFixDur(ii) = mean(curFixDur); % median gives same result
    
    xDisp = nanstd(cat(1,d.GazeX{1, d.indicesForTaskbb{ii}})); % dispersion, std
    yDisp = nanstd(cat(1,d.GazeY{1, d.indicesForTaskbb{ii}}));
    totalDispByTask(ii) = (xDisp+yDisp)/2;
end

colorsss = hsv(9);
figure(4)
for ii = 1:length(d.indicesForTaskbb)
    scatter(meanFixDur(ii),totalDispByTask(ii),'markerfacecolor',colorsss(ii,:)); hold on
end
[r,p] = corrcoef(meanFixDur,totalDispByTask);
text(.25,5.5,['r=' num2str(round(r(2),2)) ', p=' num2str(round(p(2),7))])
xlabel('Fixation duration (s)')
ylabel('Fixation dispersion (º)')
axis square
legend(d.taskList)


end