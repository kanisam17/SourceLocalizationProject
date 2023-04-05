function displayData(hPlots,subjectNameListFinal,strList,deltaPSD,freqVals,topoData,sourceData,vrbl,euDis,rangeName,refType,useMedianFlag)

%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; displaySettings.tickLengthMedium = [0.025 0];
colormap jet;
colorNames = hot(8); %colorNames([1:2,end-1:end],:) = [];
displaySettings.colorNames = colorNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titleStr = [rangeName ', ' strList{1} '(' num2str(length(subjectNameListFinal{1})) '),' strList{2} '(' num2str(length(subjectNameListFinal{2})) ')'];
cLims = [-1 1]; %cLims = [-1.25 1.25];
[~,xyz,~] = getVoxelInfo;
noseDir = '+X';
chanlocs = getMontageDetails(refType);
numGroups = length(subjectNameListFinal);
colorMrkr = {[0 0 0],[0.5 0.5 0.5]};
useMedianFlagData = true;
mtStats = cell(1, 2); % Preallocate output cell array
meuDis = cell(1,2); % Preallocate output cell array

% deltaPSD
displayAndcompareData(hPlots(1),deltaPSD,freqVals,displaySettings,cLims,1,useMedianFlag);
hold(hPlots(1),'on');
plot(hPlots(1),freqVals,zeros(1,length(freqVals)),'k');
title(hPlots(1),titleStr);

for i=1:numGroups
    % Topoplots
    subplot(hPlots(1+i));
    x = topoData{i};
    x(isnan(x)) = 999;
    topoplot_murty(x,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir',noseDir,'emarkercolors',x);
    caxis(hPlots(1+i),cLims);
    title(hPlots(1+i),strList{i});
    zoom(0.5);
    
    %%%%%%%%%%%%%%%%%%%%%% Source Localization Analysis %%%%%%%%%%%%%%%%%%%%%%%
    % Change in power in source space
    deltaPower = 10*(log10(sourceData(i).ST) - log10(sourceData(i).BL)); % Does not really work well (too noisy)
    if useMedianFlag
        mData = median(deltaPower,2);
    else
        mData = mean(deltaPower,2);
    end
    scatter3(hPlots(3+i),xyz(:,1),xyz(:,2),xyz(:,3),1,mData);
    caxis(hPlots(3+i),cLims);
    title(hPlots(3+i),['\DeltaPower ' strList{i}]);
    
end

% plot figures 
%%%%%%%% plot percentSignificant vs euclidian distance from peak. #added by kan
subplot(hPlots(6));
dataForStats = cell(1,2);
plotBinWidth = 20; plotBinLimit = 140; % across euDis dimension
nBins = plotBinLimit/plotBinWidth;
iFreq = findFreqRange(rangeName);
for numGroup = 1:2
    [h{numGroup},dataForStats{iFreq,numGroup}] = plotIndivConnData(vrbl{numGroup}(iFreq,:,:),euDis(iFreq,:),plotBinWidth,plotBinLimit,colorMrkr{numGroup},useMedianFlagData,1);
%     yline(0);
    hold on
end
if(iFreq == 1)
    set(hPlots(6), 'YDir', 'reverse');
    yticklabels([-0.1:0.2:0.5]);
    yticks([-0.1:0.2:0.5]);
    xticklabels([(plotBinWidth/2):plotBinWidth:(plotBinLimit-(plotBinWidth/2))]);
    xticks([(plotBinWidth/2):plotBinWidth:(plotBinLimit-(plotBinWidth/2))]);
else
    yticks([-0.1:0.2:0.5]);
    xticks([(plotBinWidth/2):plotBinWidth:(plotBinLimit-(plotBinWidth/2))]);
end

% statistics
% Set the confidence level
alpha = 0.05;
binEdges = 0:plotBinWidth:plotBinLimit; binned_binX = binEdges(1:end-1)+(plotBinWidth/2);
p_mat = NaN(1, nBins); % Initialize the p-value matrix
for i = iFreq  % Loop through each frequency band and bin
    for j = 1:nBins  
        data = [dataForStats{iFreq,1}(:,j); dataForStats{iFreq,2}(:,j)]; % Concatenate the data 
        group = [ones(size(dataForStats{iFreq,1},1),1); 2*ones(size(dataForStats{iFreq,2},1),1)];
        p_mat(j) = kruskalwallis(data, group, 'off'); % Perform the Kruskal-Wallis test
        if(p_mat(j) < alpha)
            baseLevel = max(nanmedian(dataForStats{iFreq,1}(:,j)),nanmedian(dataForStats{iFreq,2}(:,j)));
            scatter(binned_binX(j),baseLevel+0.12,35,'filled','p','m','HandleVisibility','off');
        end
    end
end

end
function iFreq = findFreqRange(rangeName)
if(strcmp(rangeName,'Alpha'))
    iFreq = 1;
elseif(strcmp(rangeName,'SG'))
    iFreq = 2;
elseif(strcmp(rangeName,'FG'))
    iFreq = 3;
end
end

function [h,median_binned_binY] = plotIndivConnData(mtStats,meuDis,binWidth,binLimit,colorMrkr,medianFlag,plotSwitch)
% variables list
% tStats: subjectwise t values of each voxel, subjects: 76, voxels:6239, numFreq=3, size: (3 x 6239 x 76)
% euDis: subjectwise euc. distance from max activ voxel to each voxel, dimension same as above
% binWidth: 10 as max value for euDis is 150, binned in 10 bins
% colorMrkr: used for two different groups
% plotSwitch: label to either plot or not
if ~exist('plotSwitch','var');    plotSwitch = 1;       end
connectDiscrete =  '-o';
DiscreteVisibility = 'on';
binEdges = 0:binWidth:binLimit;
nbins = length(binEdges)-1;
binned_bin_meuDis = binEdges(1:end-1)+(binWidth/2);
nsub = size(mtStats,3);
median_binned_binY = zeros(nsub,nbins); %%default mean_binY1 = zeros(nsubjects,nbins);
for i = 1:nsub
    binned_meuDis = discretize(meuDis(:,:),binEdges); %%binned_meuDis = discretize(meuDis(:,i),binEdges);
    binned_mtStats = cell(1,nbins);
    test_mtStats = mtStats(:,:,i);
    for b = 1:nbins % number of bins
        binned_mtStats{b} = test_mtStats(binned_meuDis == b);
    end
    if(medianFlag)
        median_binned_binY(i,:) = cellfun(@nanmedian,binned_mtStats); % Note:add if bin value is empty,dont assign SEM.
    else
        median_binned_binY(i,:) = cellfun(@nanmean,binned_mtStats);
    end
end
if(medianFlag)
    median_binned_binY_plot = median(median_binned_binY);
    mSEM = @(data)std(bootstrp(1000,@nanmedian,data));
else
    median_binned_binY_plot = mean(median_binned_binY);
    mSEM = @(data)std(bootstrp(1000,@nanmean,data));
end
std_binned_binY_plot = mSEM(median_binned_binY);

if(plotSwitch)
    if(length(std_binned_binY_plot)==1)
        h = errorbar(binned_bin_meuDis,median_binned_binY_plot,nan(1,8),connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
        legend('mid','old'); xlabel('eucDistance(mm)'); ylabel('deltaP')
    else
        h = errorbar(binned_bin_meuDis,median_binned_binY_plot,std_binned_binY_plot,connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
        legend('mid','old'); xlabel('eucDistance(mm)'); ylabel('deltaP')
    end
else
    h = [];
end

end

function chanlocs = getMontageDetails(refType)

capLayout = 'actiCap64';
clear cL bL chanlocs iElec electrodeList noseDir
switch refType
    case 'unipolar'
        cL = load([capLayout '.mat']);
        chanlocs = cL.chanlocs;
    case 'bipolar'
        cL = load(['bipolarChanlocs' capLayout '.mat']);
        chanlocs = cL.eloc;
end
end

function displayAndcompareData(hPlot,data,xs,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,smoothSigma,nonMatchedFlag)
if ~exist('displaySignificanceFlag','var'); displaySignificanceFlag=0;  end
if ~exist('useMedianFlag','var');           useMedianFlag=1;            end
if ~exist('smoothSigma','var');             smoothSigma=[];             end
if ~exist('nonMatchedFlag','var');          nonMatchedFlag=1;           end

if useMedianFlag
    getLoc = @(g)(squeeze(median(g,1)));
else
    getLoc = @(g)(squeeze(mean(g,1)));
end

numGroups = length(data);

if ~isempty(smoothSigma)
    windowLen = 5*smoothSigma;
    window = exp(-0.5*(((1:windowLen) - (windowLen+1)/2)/smoothSigma).^2);
    window = window/sum(window); %sqrt(sum(window.^2));
    
    for i=1:Groups
        data{i} = convn(data{i},window,'same');
    end
end

axes(hPlot);
for i=1:numGroups
    clear bootStat mData sData
    mData = getLoc(data{i}); 
    if useMedianFlag
        bootStat = bootstrp(1000,getLoc,data{i});
        sData = std(bootStat);
    else
        sData = std(data{i},[],1)/sqrt(size(data{i},1));
    end
    
    patch([xs';flipud(xs')],[mData'-sData';flipud(mData'+sData')],displaySettings.colorNames(i,:),'linestyle','none','FaceAlpha',0.4);
    hold on;
    plot(xs,mData,'color',displaySettings.colorNames(i,:),'linewidth',1);
end

set(gca,'fontsize',displaySettings.fontSizeLarge);
set(gca,'TickDir','out','TickLength',displaySettings.tickLengthMedium);

if exist('yLims','var') && ~isempty(yLims)
    ylim(yLims);
else
    yLims = ylim;
end

if displaySignificanceFlag % Do significance Testing
    
    allData = [];
    allIDs = [];
    for j=1:numGroups
        allData = cat(1,allData,data{j});
        allIDs = cat(1,allIDs,j+zeros(size(data{j},1),1));
    end
       
   for i=1:length(xs)
       if useMedianFlag
           p=kruskalwallis(allData(:,i),allIDs,'off');
       else
           if nonMatchedFlag
               [~,p]=ttest2(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
           else
               [~,p]=ttest(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
           end
       end
       % Get patch coordinates
       yVals = yLims(1)+[0 0 diff(yLims)/20 diff(yLims)/20];
       
       clear xMidPos xBegPos xEndPos
       xMidPos = xs(i);
       if i==1
           xBegPos = xMidPos;
       else
           xBegPos = xMidPos-(xs(i)-xs(i-1))/2; 
       end
       if i==length(xs)
           xEndPos = xMidPos; 
       else
           xEndPos = xMidPos+(xs(i+1)-xs(i))/2; 
       end
       clear xVals; xVals = [xBegPos xEndPos xEndPos xBegPos]';
       
       if (p<0.05)
           patch(xVals,yVals,'k','linestyle','none');
       end
       if (p<0.01)
           patch(xVals,yVals,'g','linestyle','none');
       end
   end
end
end