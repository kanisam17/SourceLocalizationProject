clear; clc;
displaySingleGroupFlag=0; % When a single group is provided, stimulus and baseline epochs are compared
useMedianFlag=1;

thres=-10; % only subjects who have delta power above this (in dB) are selected. Set to a low value such as -inf to take all subjects
useCommonSubjectsFlag=1; % if set to 1, only subjects for which delta power is more than threshold for all frequencies are chosen

%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; displaySettings.tickLengthMedium = [0.025 0];
% colormap magma;
colormap jet;
colorNames = hot(8); colorNames([1:3,end-2:end],:) = [];
displaySettings.colorNames = colorNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mandatory fixed options
%folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];

% Choose one of these options
refType = 'unipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjects = getGoodSubjectsProjectwise(projectName,1);
uniqueSubjectNames0 = getGoodFileNamesForSubjects(goodSubjects{1});

%%%%%%%%%%%%%% Find indices for which the correct capType was used %%%%%%%%
capTypeToUse = 'actiCap64';
goodIndices = [];
for i=1:length(uniqueSubjectNames0)
    [expDates,~,capType,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,uniqueSubjectNames0{i},protocolType);
    if usableDataFlag && ~isempty(expDates) && strcmp(capType{1},capTypeToUse)
        goodIndices = cat(2,goodIndices,i);
    end
end
disp([num2str(length(goodIndices)) ' subjects with correct capType chosen for further analysis']);
uniqueSubjectNames = uniqueSubjectNames0(goodIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%% Load Power Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma1Range = [20 34]; gamma2Range = [36 66]; alphaRange = [8 12];
spatialFrequenciesToRemove=[];
dataForDisplay = combineAnalyzedData(pwd,uniqueSubjectNames,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,0);

%%%%%%%%%%%%%%%%%%%%%%%% Find Useful Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);
healthyPos = strcmp(cdrList,'HV');
ageGroup1Pos = (ageList<65) & healthyPos; strList{1} = 'mid';
ageGroup2Pos = (ageList>=65) & healthyPos; strList{2} = 'old';
    
numFreqRanges = length(dataForDisplay.rangeNames);

hPlots = getPlotHandles(numFreqRanges,3,[0.05 0.05 0.4 0.9],0.02,0.05);

dataDeltaPSD = 10*(dataForDisplay.logSTPowerVsFreqAllSubjects - dataForDisplay.logBLPowerVsFreqAllSubjects);
cLims = [-1.25 1.25];
noseDir = '+X';
chanlocs = getMontageDetails(refType);

% Select common good subjects if needed
goodSubjectPosAll = zeros(size(dataForDisplay.powerDBAllSubjects));
for i=1:numFreqRanges
    if strcmp(dataForDisplay.rangeNames{i},'Alpha')
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)<-thres;
    else
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)>thres;
    end
end
if useCommonSubjectsFlag 
    goodSubjectPosCommon = all(goodSubjectPosAll,2);
end

plotPosition = [2 3 1]; % The data are saved as SG, FG, Alpha. But we want to plot such in order: Alpha, SG, FG.

for i=1:numFreqRanges
    if useCommonSubjectsFlag
        goodSubjectPos = goodSubjectPosCommon;
    else
        goodSubjectPos = goodSubjectPosAll(:,i);
    end
    
    clear goodPos subjectNameListFinal
    gp1 = ageGroup1Pos & goodSubjectPos';
    gp2 = ageGroup2Pos & goodSubjectPos';
    
    if displaySingleGroupFlag
        goodPos{1} = gp1 | gp2;
        subjectNameListFinal{1} = uniqueSubjectNames(goodPos{1});
        deltaPSD{1} = dataDeltaPSD(goodPos{1},:);
        titleStr = ['Combined(' num2str(length(subjectNameListFinal{1})) ')'];
        displaySignificanceFlag=0;
    else
        goodPos{1} = gp1; %#ok<*UNRCH>
        goodPos{2} = gp2;
        subjectNameListFinal{1} = uniqueSubjectNames(gp1);
        subjectNameListFinal{2} = uniqueSubjectNames(gp2);
        deltaPSD{1} = dataDeltaPSD(goodPos{1},:);
        deltaPSD{2} = dataDeltaPSD(goodPos{2},:);
        titleStr = [strList{1} '(' num2str(length(subjectNameListFinal{1})) '),' strList{2} '(' num2str(length(subjectNameListFinal{2})) ')'];
        displaySignificanceFlag=1;
    end

    plotPos = plotPosition(i);
    displayAndcompareData(hPlots(plotPos,1),deltaPSD,dataForDisplay.freqVals,displaySettings,cLims,displaySignificanceFlag,useMedianFlag);
    hold(hPlots(plotPos,1),'on');
    plot(hPlots(plotPos,1),dataForDisplay.freqVals,zeros(1,length(dataForDisplay.freqVals)),'k');
    title(hPlots(plotPos,1),[dataForDisplay.rangeNames{i} ', ' titleStr]);
    
    numGroups = length(subjectNameListFinal);   
    for j=1:numGroups        
        axes(hPlots(plotPos,1+j)); %#ok<LAXES>
        if useMedianFlag
            x = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,goodPos{j}),3));
        else
            x = squeeze(nanmean(dataForDisplay.powerDBTopoAllSubjects(:,i,goodPos{j}),3));
        end
        x(isnan(x)) = 999;
        topoplot_murty(x,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir',noseDir,'emarkercolors',x);
        caxis(cLims);
    end
end

%%%%%%%%%%%%%%%%%%%%%% Source Localization Analysis %%%%%%%%%%%%%%%%%%%%%%%
disp('Getting sLORETA data');
folderLORETA = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\Kanishka_SourceLocalizationProject\data\sLORETA_Thres10';
[allDataBL,allDataST] = getLORETAData(subjectNameListFinal,strList,folderLORETA);
[posList,xyz,areaList] = getVoxelInfo;
numAreas = length(areaList);
colorNamesAreas = jet(numAreas);

hPlotsSource = getPlotHandles(numFreqRanges,5,[0.5 0.05 0.475 0.9],0.02,0.05);
for i=1:numFreqRanges
    
    % First, we compare stim versus baseline responses, separately for each group
    for j=1:2
        compareData{1} = squeeze(allDataBL{j}(:,i,:));
        compareData{2} = squeeze(allDataST{j}(:,i,:));

        deltaPower = 10*(log10(compareData{2}) - log10(compareData{1}));
        if useMedianFlag
            mData = median(deltaPower);
            %sData = getSEMedian(deltaPower);
        else
            mData = mean(deltaPower);
            sData = std(deltaPower);
        end
        scatter3(hPlotsSource(i,2*(j-1)+1),xyz(:,1),xyz(:,2),xyz(:,3),1,mData);
        caxis(hPlotsSource(i,2*(j-1)+1),cLims);
        
        startPos=0;
        hold(hPlotsSource(i,2*j),'on');
        for k=1:numAreas
            mDataTMP = (mData(posList{k}));
            numEntries = length(mDataTMP);
            plot(hPlotsSource(i,2*j),startPos+(1:numEntries),mDataTMP,'color',colorNamesAreas(k,:));
            startPos = startPos+numEntries;
        end
        ylim(hPlotsSource(i,2*j),cLims);
    end
    % displayAndcompareData(hPlotsSource(i,3),compareData,1:6239,displaySettings,[],displaySignificanceFlag,useMedianFlag);
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
    
    for i=1:numGroups
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