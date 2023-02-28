% ToDo
% Display
% 0. Display - show better brain maps. Currently using scatter3.
% 1. Display - 8th column - plot errorbars properly as swarm plots.
% 2. Display - the 9th column show show fraction as a function from the
% voxel that has maximum activation. 

% Analysis
% 1. Subject power matching
% 2. Find electrodes that are bad in many subjects and remove those common
% bad electrodes. Do LORETA analysis on the subset of subjects who have remaining
% good electrodes
% 3. Number of trials matching - since stats may depend on the number of
% trials, control for that by having the same number of trials for all subjects.
% 4. If there are enough subjects for Case/Control - do that comparison.

clear; clc;

thres=-10; % only subjects who have delta power above this (in dB) are selected. Set to a low value such as -inf to take all subjects
useCommonSubjectsFlag=1; % if set to 1, only subjects for which delta power is more than threshold for all frequencies are chosen
useMedianFlag=1;
% folderSourceString = pwd; % input for powermatching
% folderLORETA = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\randomTrails\sLORETA_Thres10';
folderLORETA = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\sourceData\LORETA\data\Age'; % Folder where the output of LORETA is saved;

% powerMatchedSubjectList = load ('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\matchedSubjectNameList_FG_new.mat');
powerMatchedSubjectList = load ('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\matchedSubjectNameList_A.mat');

caseList = load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\ADGammaProjectCodes\caseAgeMatchedSubjectList.mat');
% folderSourceString = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10';
xyz = xlsread ('voxelInfo.xlsx');

% powerMatchedSubjectList.matchedSubjectNameLists = powerMatchedSubjectList.powerMatchedSubjectNameLists;
powerMatchedSubjectList.matchedSubjectNameLists = powerMatchedSubjectList.matchedSubjectNameList_A; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mandatory fixed options
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
goodIndices = []; protocolNamesUnique = {};expDatesUnique = {};
for i=1:length(uniqueSubjectNames0)
    [expDates,protocolNames,capType,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,uniqueSubjectNames0{i},protocolType);
    if usableDataFlag && ~isempty(expDates) && strcmp(capType{1},capTypeToUse)
        goodIndices = cat(2,goodIndices,i); 
    end
    protocolNames0 = {protocolNames};
    try
    protocolNamesUnique = cat(2,protocolNamesUnique,protocolNames0);
    catch
    end
    expDate0 = {expDates};
    try
    expDatesUnique = cat(2,expDatesUnique,expDate0);
    catch
    end
end

for i = 1:size(protocolNamesUnique,2)
    protocolNamesUnique{i} = protocolNamesUnique{i}';
end


disp([num2str(length(goodIndices)) 'subjects with correct capType chosen for further analysis']);
uniqueSubjectNames = uniqueSubjectNames0(goodIndices);


%%%%%%%%%%%%%%%%%%%%%%%%%% Creating random trials index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analyzedDataLocation = 'C:\Users\Liza\Documents\TLSAEEGProject\SourceLocalizationProject\23Nov2022\SourceLocalizationProject\analyzedDataRandTrials';
analyzedDataFolder = fullfile(analyzedDataLocation, projectName, protocolType);
[trialIdxAllSub0,~,subNameIdx,~] = getTrialIdxAndGoodSubject64Elecs(analyzedDataLocation,projectName,protocolType,uniqueSubjectNames);
trialIdxAllSub = trialIdxAllSub0(subNameIdx); %trialIdxAllSub is index of random trials selected for uniqueSubjectNames.



%%%%%%%%%%%%%%%%%%%%%%%%%% Load Power Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma1Range = [20 34]; gamma2Range = [36 66]; alphaRange = [8 12];
spatialFrequenciesToRemove=[];
%dataForDisplay = combineAnalyzedData(pwd,uniqueSubjectNames,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,0);
dataForDisplay = combineAnalyzedData(pwd,uniqueSubjectNames,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,0);

% Here data is saved in the order [SG, FG, alpha]. Change to [alpha SG FG];
newOrderList = [3 1 2];
dataForDisplay.rangeNames = dataForDisplay.rangeNames(newOrderList);
dataForDisplay.powerDBAllSubjects = dataForDisplay.powerDBAllSubjects(:,newOrderList);
dataForDisplay.powerDBTopoAllSubjects = dataForDisplay.powerDBTopoAllSubjects(:,newOrderList,:);

%%%%%%%%%%%%%%%%%%%%%%%% Find Useful Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);
healthyPos = strcmp(cdrList,'HV');
malePos = strcmp(genderList, 'M'); femalePos = strcmp(genderList, 'F');
%%
conditions = {'PowdB_10','trialsControl','powerMatch'};
%% mid or old
% ageGroup1Pos = (ageList<65) & healthyPos;
% ageGroup2Pos = (ageList>=65) & healthyPos;
ageGroup1Pos = (ageList<65);
ageGroup2Pos = (ageList>=65);
%% Case or control %section Added by kan
casePos1 = strcmp(cdrList,'MCI'); 
casePos2 = strcmp(cdrList,'AD');
casePos = casePos1 | casePos2;

controlList = cell(1,length(caseList.subjectNameListMatched{1,1}));
for i = 1:length(caseList.subjectNameListMatched{1,1})
    randomIndex = randi(length(caseList.subjectNameListMatched{1,1}{i}));
    controlList{i} = caseList.subjectNameListMatched{1,1}{i}{randomIndex};
end
[controlPos,controlidx] = ismember(uniqueSubjectNames,controlList);
%%
%% trials condition trials>180
% ageGroup1idx = find((ageList<65) & healthyPos); strList{1} = 'mid';
% ageGroup2idx = find((ageList>=65) & healthyPos); strList{2} = 'old';
ageGroup1idx = find(ageGroup1Pos); strList{1} = 'mid'; %strList{1} = 'control'; %
ageGroup2idx = find(ageGroup2Pos); strList{2} = 'old'; %strList{2} = 'case'; %
listNumStimuli = trialsAvailableEachSubject();
ageGroup1idx = ageGroup1idx(listNumStimuli{1}>180);
ageGroup2idx = ageGroup2idx(listNumStimuli{2}>180);
trialMatchPos = false(1,237);
trialMatchPos(ageGroup1idx) = 1;
trialMatchPos(ageGroup2idx) = 1;
%
% To find powermatched mid and old list
[powerMatch1Pos,~] =ismember(uniqueSubjectNames, powerMatchedSubjectList.matchedSubjectNameLists{1});strList{1} = 'mid';
[powerMatch2Pos,~] =ismember(uniqueSubjectNames, powerMatchedSubjectList.matchedSubjectNameLists{2});strList{2} = 'old';
powerMatch = powerMatch1Pos | powerMatch2Pos;

%% 
numFreqRanges = length(dataForDisplay.rangeNames);
dataDeltaPSD = 10*(dataForDisplay.logSTPowerVsFreqAllSubjects - dataForDisplay.logBLPowerVsFreqAllSubjects);

% Select common good subjects if needed
goodSubjectPosAll = zeros(size(dataForDisplay.powerDBAllSubjects));
for i=1:numFreqRanges
    if strcmp(dataForDisplay.rangeNames{i},'Alpha') % Alpha being lowest powerDB selected to cover all subjects.
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)<-thres;
    else
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)>thres;
    end
end
if useCommonSubjectsFlag 
    goodSubjectPosCommon = all(goodSubjectPosAll,2);
end

%%
% Generate plots with 9 columns
% 1) deltaPSD, 2) topoBL, 3) topoST, 
% 4-5) sourceBL & sourceST - change in power
% 6-7) sourceBL & sourceST - tStat
% 8) percentSignificant - barplot
% 9) percentSignificant vs distance from peak

hPlots = getPlotHandles(numFreqRanges,9,[0.05 0.05 0.9 0.9],0.02,0.05);

for i=1:numFreqRanges
    if useCommonSubjectsFlag
        goodSubjectPos = goodSubjectPosCommon;
    else
        goodSubjectPos = goodSubjectPosAll(:,i); %#ok<*UNRCH>
    end
    
    clear goodPos subjectNameListFinal
    %     gp1 = ageGroup1Pos & goodSubjectPos';
    %     gp2 = ageGroup2Pos & goodSubjectPos';
    %     gp1 = ageGroup1Pos & goodSubjectPos'& healthyPos & malePos;% & powerMatch;% & trialMatchPos;
    %     gp2 = ageGroup2Pos & goodSubjectPos'& healthyPos & malePos;% powerMatch;% & trialMatchPos;
    gp1 = ageGroup1Pos & goodSubjectPos' & healthyPos;% & powerMatch;% & trialMatchPos;
    gp2 = ageGroup2Pos & goodSubjectPos' & healthyPos;% & powerMatch;% & trialMatchPos;
    
    subjectNameListFinal{1} = uniqueSubjectNames(gp1);
    subjectNameListFinal{2} = uniqueSubjectNames(gp2);
    trialIdxListFinal{1} = trialIdxAllSub(gp1);% Added by Kan
    trialIdxListFinal{2} = trialIdxAllSub(gp2);% Added by Kan
    protocolNamesFinal{1} = protocolNamesUnique(gp1);% Added by Kan
    protocolNamesFinal{2} = protocolNamesUnique(gp2);% Added by Kan
    expDatesFinal{1} = expDatesUnique(gp1);% Added by Kan
    expDatesFinal{2} = expDatesUnique(gp2);% Added by Kan
%     goodProtFlagFinal{1} = goodProtFlag(gp1);% Added by Kan
%     goodProtFlagFinal{2} = goodProtFlag(gp2);% Added by Kan
    
    deltaPSD{1} = dataDeltaPSD(gp1,:);
    deltaPSD{2} = dataDeltaPSD(gp2,:);
    
    
    if useMedianFlag
        topoData{1} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp1),3));
        topoData{2} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp2),3));
    else
        topoData{1} = squeeze(nanmean(dataForDisplay.powerDBTopoAllSubjects(:,i,gp1),3));
        topoData{2} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp2),3));
    end
    
    [allDataBL,allDataST,alltStats,allpVals] = getLORETAData(subjectNameListFinal,strList,folderLORETA);
    for j=1:2
        sourceData(j).BL = squeeze(allDataBL{j}(:,i,:)); %#ok<*SAGROW>
        sourceData(j).ST = squeeze(allDataST{j}(:,i,:));
%         sourceData(j).DeltaP = squeeze(allDataDeltaP{j}(:,i,:));
        sourceData(j).tStats = squeeze(alltStats{j}(:,i,:));
        sourceData(j).pVals = squeeze(allpVals{j}(:,i,:));
    end
    [euDis,tStats,dataDeltaP] = displayData(hPlots(i,:),subjectNameListFinal,strList,deltaPSD,dataForDisplay.freqVals,topoData,sourceData,dataForDisplay.rangeNames{i},refType,useMedianFlag,folderLORETA,xyz);
end


%%  To find powermatched mid and old list
% powerMatchedSubjectNameLists = getPowerMatchedSubjectList(folderSourceString,subjectNameListFinal,projectName,refType,protocolType);