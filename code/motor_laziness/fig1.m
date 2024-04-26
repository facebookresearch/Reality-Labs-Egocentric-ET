function [] = fig1(d)
% Copyright (c) Meta Platforms, Inc. and affiliates.
% fig1

d = d.GazeInHead;

%% DATA: boostrap to get equal data across tasks

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
xBins = linspace(-40,40,60);
yBins = linspace(-25,25,60);

% generate fixation position and duration maps
for ii = 1:length(xBins)-1
    for jj = 1:length(yBins)-1
        justDataInBin = fixAllll(posAlllll(:,1)>=xBins(ii) & posAlllll(:,1)<xBins(ii+1) & posAlllll(:,2)>=yBins(jj) & posAlllll(:,2)<yBins(jj+1));
        meanmap1(jj,ii) = nanmean(justDataInBin);
        semMap(jj,ii) = nanstd(justDataInBin)./sqrt(length(justDataInBin));
        nPointsmap1(jj,ii) = length(justDataInBin);
    end
end


figure(1);
subplot(3,1,2);
imagesc(xBins,yBins,nPointsmap1)
axis equal; hline(0,'-k'); vline(0,'-k');
set(gca,'YDir','normal')
colormap(flipud(parula))
xlim([-40 40])
ylim([-25 25])
colorbar

subplot(3,1,3);
imagesc(xBins,yBins,meanmap1,'AlphaData',nPointsmap1>=22 ) % 22 obeys the 2.2 data:bootstrappedData ratio, equivalent to removing 10 datapoints
axis equal; hline(0,'-k'); vline(0,'-k');
set(gca,'YDir','normal')
colormap(flipud(parula))
caxis([.1 .8])
ylim([-25 25])
xlim([-40 40])
colorbar


% make cartoon
subplot(3,1,1); 

x1 = linspace(-45,45,60*3); x2=linspace(-25,25,60*3);
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];
mu = [0,0];
Sigma1 = eye(2).*110;
Sigma2 = eye(2).*50;
y = mvnpdf(X,mu,Sigma1) + 0.2.*mvnpdf(X,mu,Sigma2);
y = reshape(y,length(x2),length(x1));
imagesc(x1,x2,y);
 hline(0,'-k'); vline(0,'-k');
set(gca,'YDir','normal')
axis equal
colormap(flipud(parula))
caxis([0 max(max(y))])
ylim([-25 25])
xlim([-45 45])



