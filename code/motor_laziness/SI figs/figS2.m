function [] = figS2(d)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% figS2

d = d.GazeInHead;


%% spatial binning
xBins = linspace(-40,40,60);
yBins = linspace(-25,25,60);

% per task maps
figure;
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

    subplot(3,3,bb);
    imagesc(xBins,yBins,meanmapPerTask{bb},'AlphaData',nPointsmapPerTask{bb}>=1 ); % uncomment to plot fixation duration map
    %imagesc(xBins,yBins,nPointsmapPerTask{bb});%'AlphaData',nPointsmapPerTask{bb}./max(max(nPointsmapPerTask{bb}))); % uncomment to plot fixation location map
    set(gca,'YDir','normal')
    %clim([0 1.5]) % same color limit for each plot
    title(d.taskList{bb})
    %axis equal
    hline(0,'-k'); vline(0,'-k');
    colormap(flipud(parula))
    caxis([.1 .8])
    %ylim([-25 25])
    xlim([-40,40])
    if bb == 7
        xlabel('Horizontal position (º)')
        ylabel('Vertical position (º)')
    end


end


end