function [euDis,tStats,dataDeltaP,mData] = displayData(hPlots,subjectNameListFinal,strList,deltaPSD,freqVals,topoData,sourceData,rangeName,refType,useMedianFlag,folderLORETA,xyz,idxFPTO)

%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; displaySettings.tickLengthMedium = [0.025 0];
% colormap magma;
colormap jet;
colorNames = hot(8); %colorNames([1:2,end-1:end],:) = [];
displaySettings.colorNames = colorNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titleStr = [rangeName ', ' strList{1} '(' num2str(length(subjectNameListFinal{1})) '),' strList{2} '(' num2str(length(subjectNameListFinal{2})) ')'];
cLims = [-1 1]; %cLims = [-1.25 1.25];

% deltaPSD
displayAndcompareData(hPlots(1),deltaPSD,freqVals,displaySettings,cLims,1,useMedianFlag);
hold(hPlots(1),'on');
plot(hPlots(1),freqVals,zeros(1,length(freqVals)),'k');
title(hPlots(1),titleStr);

% Topoplots
noseDir = '+X';
chanlocs = getMontageDetails(refType);
numGroups = length(subjectNameListFinal);
for i=1:numGroups
    axes(hPlots(1+i)); %#ok<LAXES>
    x = topoData{i};
    x(isnan(x)) = 999;
    topoplot_murty(x,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir',noseDir,'emarkercolors',x);
    caxis(cLims);
    title(strList{i});
end

%%%%%%%%%%%%%%%%%%%%%% Source Localization Analysis %%%%%%%%%%%%%%%%%%%%%%%
[posList,xyz,areaList] = getVoxelInfo;

% Change in power in source space
for i=1:numGroups
    deltaPower = 10*(log10(sourceData(i).ST) - log10(sourceData(i).BL)); % Does not really work well (too noisy)
    if useMedianFlag
        mData = median(deltaPower);
    else
        mData = mean(deltaPower);
    end
    scatter3(hPlots(3+i),xyz(:,1),xyz(:,2),xyz(:,3),1,mData);
    caxis(hPlots(3+i),cLims);
    title(hPlots(3+i),['\DeltaPower ' strList{i}]);
end

% Plot t-stats
for i=1:numGroups
    statVals = sourceData(i).tStats;
    if useMedianFlag
        mData = median(statVals);
    else
        mData = mean(statVals);
    end
    scatter3(hPlots(5+i),xyz(:,1),xyz(:,2),xyz(:,3),1,mData); %%
    caxis(hPlots(5+i),cLims);
    title(hPlots(5+i),['tStats ' strList{i}]);
end

% plot fraction of significant voxels for each area
numAreas = length(areaList);
pThreshold = 0.05;

allFractionLists = cell(1,numGroups);
for i=1:numGroups
    allpVals = sourceData(i).pVals;
    numSubjects = size(allpVals,1);
        
    fractionList = zeros(numSubjects,numAreas);
    for s=1:numSubjects
        pVals = squeeze(allpVals(s,:));
        for k=1:numAreas
            pDataTMP = (pVals(posList{k}));
            fractionList(s,k) = length(find(pDataTMP<pThreshold))/length(pDataTMP);
        end
    end
    allFractionLists{i} = fractionList;
end

del=0.3; %del=0.3;
hold(hPlots(8),'on');

pList = zeros(1,numAreas);
areaStr = cell(1,numAreas);
for i=1:numAreas
    areaStr{i} = areaList{i}(1);
    x1 = allFractionLists{1}(:,i);
    x2 = allFractionLists{2}(:,i);
    z1 = zscore(allFractionLists{1}(:,i));
    z2 = zscore(allFractionLists{2}(:,i));
    
    
    if useMedianFlag
        mX1 = median(x1); sX1 = getSEMedian(x1);
        mX2 = median(x2); sX2 = getSEMedian(x2);
        pList(i) = ranksum(x1,x2);
    else
        mX1 = mean(x1); sX1 = std(x1)/sqrt(length(x1));
        mX2 = mean(x2); sX2 = std(x2)/sqrt(length(x2));
        [~,pList(i)] = ttest2(x1,x2);
    end
    %Log scale plotting
    %x1 = 10.log10(x1); x2 = 10.log10(x2); mX1 = 10.log10(mX1); mX2 = 10.log10(mX2);
    
    bar(i-del/2,mX1,'FaceAlpha',0.8,'facecolor',displaySettings.colorNames(2,:),'barwidth',2*del/3,'Parent',hPlots(8));
    bar(i+del/2,mX2,'FaceAlpha',0.8,'facecolor',displaySettings.colorNames(3,:),'barwidth',2*del/3,'Parent',hPlots(8));
    
    
    swarmchart(ones(1,length(x1))*i-del/2,x1,11,displaySettings.colorNames(1,:),'filled','MarkerFaceAlpha',0.3,'Parent',hPlots(8)); %Added by Kan
    swarmchart(ones(1,length(x2))*i+del/2,x2,11,displaySettings.colorNames(4,:),'filled','MarkerFaceAlpha',0.3,'Parent',hPlots(8)); %Added by Kan
   
    errorbar(hPlots(8),i-del/2,mX1,sX1,'color',displaySettings.colorNames(1,:));
    errorbar(hPlots(8),i+del/2,mX2,sX2,'color',displaySettings.colorNames(4,:));
    
    title(hPlots(8),['Activated/Total voxel in lobes']);
end
set(hPlots(8),'XTick',1:numAreas,'XTickLabel',areaStr);
disp([rangeName ', p=' num2str(pList)]);


% % plot percentSignificant vs euclidian distance from peak. #added by kan
% 
% [euDis, tStats, meuDis, dataDeltaP, mtStats]= euDistanceForLORETA(folderLORETA,subjectNameListFinal,strList,xyz,1);
[euDis,tStats,dataDeltaP]= euDistanceForLORETA(folderLORETA,subjectNameListFinal,strList,xyz,0,idxFPTO);

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