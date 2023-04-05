% Analysis
% 1. Subject power matching
% 2. Find electrodes that are bad in many subjects and remove those common
% bad electrodes. Do LORETA analysis on the subset of subjects who have remaining
% good electrodes
% 3. Number of trials matching - since stats may depend on the number of
% trials, control for that by having the same number of trials for all subjects.
% 4. If there are enough subjects for Case/Control - do that comparison.

clear; clc;

deltaThresSubjSelect=-10; % only subjects who have delta power above this (in dB) are selected. Set to a low value such as -inf to take all subjects
useCommonSubjectsFlag=0; % if set to 1, only subjects for which delta power is more than threshold for all frequencies are chosen
useMedianFlag=1;
% folderSourceString = pwd; % input for powermatching
% folderLORETA = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\randomTrails\sLORETA_Thres10';
% folderLORETA = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\sourceData\LORETA\data\Age'; % Folder where the output of LORETA is saved;
% InterpolatedData
% folderLORETA = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10\interpolatedData\decimatedData\sourceData\LORETA\data\text';
folderLORETA = 'Z:\Projects\Kanishka_SourceLocalizationProject\data\interpolatedData';
% powerMatchedSubjectList = load ('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\matchedSubjectNameList_FG_new.mat');
% powerMatchedSubjectList = load ('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\matchedSubjectNameList_SG.mat');
% powerMatchedSubjectList.matchedSubjectNameLists = powerMatchedSubjectList.powerMatchedSubjectNameLists;

% caseList = load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\ADGammaProjectCodes\caseAgeMatchedSubjectList.mat');

xyz = xlsread ('voxelInfo.xlsx');
lobes = readtable('voxelInfo');
lobesName = table2array(lobes(:,4));
idxFPTO = ismember(lobesName, {'Limbic Lobe', 'Sub-lobar', 'Temporal Lobe'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];
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
    try
        protocolNamesUnique = cat(2,protocolNamesUnique,{protocolNames});
        expDatesUnique = cat(2,expDatesUnique,{expDates});
    catch
    end
end

disp([num2str(length(goodIndices)) 'subjects with correct capType chosen for further analysis']);
uniqueSubjectNames = uniqueSubjectNames0(goodIndices);


%%%%%%%%%%%%%%%%%%%%%%%%%% Creating random trials index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyzedDataLocation = 'C:\Users\Liza\Documents\TLSAEEGProject\SourceLocalizationProject\23Nov2022\SourceLocalizationProject\analyzedDataRandTrials';
% analyzedDataFolder = fullfile(analyzedDataLocation, projectName, protocolType);
% [trialIdxAllSub0,~,subNameIdx,~] = getTrialIdxAndGoodSubject64Elecs(analyzedDataLocation,projectName,protocolType,uniqueSubjectNames);
% trialIdxAllSub = trialIdxAllSub0(subNameIdx); %trialIdxAllSub is index of random trials selected for uniqueSubjectNames.
prevAnalyzedDataLoc = 'Z:\Projects\Kanishka_SourceLocalizationProject\data';
%%%%%%%%%%%%%%%%%%%%%%%%%% Load Power Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma1Range = [20 34]; gamma2Range = [36 66]; alphaRange = [8 12]; spatialFrequenciesToRemove=[];
dataForDisplay = combineAnalyzedData(prevAnalyzedDataLoc,uniqueSubjectNames,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,0);
% Here data is saved in the order [SG, FG, alpha]. Change to [alpha SG FG];
newOrderList = [3 1 2];dataForDisplay.rangeNames = dataForDisplay.rangeNames(newOrderList);
dataForDisplay.powerDBAllSubjects = dataForDisplay.powerDBAllSubjects(:,newOrderList);
dataForDisplay.powerDBTopoAllSubjects = dataForDisplay.powerDBTopoAllSubjects(:,newOrderList,:);

%%%%%%%%%%%%%%%%%%%%%%%% Find Useful Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);
healthyPos = strcmp(cdrList,'HV');
malePos = strcmp(genderList, 'M'); femalePos = strcmp(genderList, 'F');

% mid or old
ageGroup1Pos = (ageList<65) & healthyPos;
ageGroup2Pos = (ageList>=65) & healthyPos;
% Case or control %section Added by kan
casePos1 = strcmp(cdrList,'MCI'); 
casePos2 = strcmp(cdrList,'AD');
casePos = casePos1 | casePos2;

comparisonCond = 'midvsold';

if(strcmp(comparisonCond,'midvsold'))
conditions = {'PowdB_10','trialsControl','powerMatch'};
pCondition = 1; % [1 2]; [1 2 3];
strList{1} = 'mid'; strList{2} = 'old';
%% trials condition trials>180
if(strcmp(conditions{pCondition},'trialsControl'))
ageGroup1idx = find(ageGroup1Pos); strList{1} = 'mid'; %strList{1} = 'control'; %
ageGroup2idx = find(ageGroup2Pos); strList{2} = 'old'; %strList{2} = 'case'; %
listNumStimuli = trialsAvailableEachSubject();
ageGroup1idx = ageGroup1idx(listNumStimuli{1}>180);
ageGroup2idx = ageGroup2idx(listNumStimuli{2}>180);
trialMatchPos = false(1,length(uniqueSubjectNames));
trialMatchPos(ageGroup1idx) = 1; trialMatchPos(ageGroup2idx) = 1;
end
%% To find powermatched mid and old list
% if(strcmp(conditions{pCondition},'powerMatch'))
% [powerMatch1Pos,~] =ismember(uniqueSubjectNames, powerMatchedSubjectList.matchedSubjectNameLists{1});strList{1} = 'mid';
% [powerMatch2Pos,~] =ismember(uniqueSubjectNames, powerMatchedSubjectList.matchedSubjectNameLists{2});strList{2} = 'old';
% powerMatch = powerMatch1Pos | powerMatch2Pos;
% end

else % 'casevscontrol'
    strList{1} = 'control'; strList{2} = 'case';
end


%%
numFreqRanges = length(dataForDisplay.rangeNames);
dataDeltaPSD = 10*(dataForDisplay.logSTPowerVsFreqAllSubjects - dataForDisplay.logBLPowerVsFreqAllSubjects);
goodSubjectPosAll = zeros(size(dataForDisplay.powerDBAllSubjects));
for i=1:numFreqRanges
    if strcmp(dataForDisplay.rangeNames{i},'Alpha') % Alpha being lowest powerDB selected to cover all subjects.
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)<-deltaThresSubjSelect;
    else
        goodSubjectPosAll(:,i) = dataForDisplay.powerDBAllSubjects(:,i)>deltaThresSubjSelect;
    end
end
goodSubjectPosCommon = all(goodSubjectPosAll,2);


% Generate plots
hPlots = getPlotHandles(numFreqRanges,6,[0.05 0.05 0.9 0.9],0.02,0.05);
for i=1:numFreqRanges
    if useCommonSubjectsFlag
        goodSubjectPos = goodSubjectPosCommon;
    else
        goodSubjectPos = goodSubjectPosAll(:,i); %#ok<*UNRCH>
    end
    if(strcmp(comparisonCond,'midvsold'))
        clear gp1 gp2 subjectNameListFinal
        gp1 = ageGroup1Pos & goodSubjectPos' & healthyPos;
        gp2 = ageGroup2Pos & goodSubjectPos' & healthyPos;
        if(strcmp(conditions{pCondition},'powerMatch'))
            gp1 = gp1 & powerMatch; gp2 = gp2 & powerMatch;
        end
        if(strcmp(conditions{pCondition},'trialsControl'))
            gp1 = gp1 & trialMatchPos; gp2 = gp2 & trialMatchPos;
        end
        subjectNameListFinal{1} = uniqueSubjectNames(gp1);
        subjectNameListFinal{2} = uniqueSubjectNames(gp2);
        
        deltaPSD{1} = dataDeltaPSD(gp1,:);
    deltaPSD{2} = dataDeltaPSD(gp2,:);
    if useMedianFlag
        topoData{1} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp1),3));
        topoData{2} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp2),3));
    else
        topoData{1} = squeeze(nanmean(dataForDisplay.powerDBTopoAllSubjects(:,i,gp1),3));
        topoData{2} = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,i,gp2),3));
    end
    
    else % 'casevscontrol'
        subjectNameListFinal = getCaseControlList(uniqueSubjectNames,cdrList,ageList,genderList);
        
%         deltaPSD
        
%         topoData
    end
%     trialIdxListFinal{1} = trialIdxAllSub(gp1);% Added by Kan
%     trialIdxListFinal{2} = trialIdxAllSub(gp2);% Added by Kan
    
    
    [allDataBL,allDataST,allDataDeltaP,alltStats,allpVals] = getLORETAData(subjectNameListFinal,strList,folderLORETA);
    vrbl = allDataDeltaP;
    euDis = euDistanceForLORETA(vrbl,xyz,useMedianFlag);
    for j=1:2
        sourceData(j).BL = squeeze(allDataBL{j}(i,:,:)); %#ok<*SAGROW>
        sourceData(j).ST = squeeze(allDataST{j}(i,:,:));
        sourceData(j).DeltaP = squeeze(allDataDeltaP{j}(i,:,:));
        sourceData(j).tStats = squeeze(alltStats{j}(i,:,:));
        sourceData(j).pVals = squeeze(allpVals{j}(i,:,:));
    end
    displayData(hPlots(i,:),subjectNameListFinal,strList,deltaPSD,dataForDisplay.freqVals,topoData,sourceData,vrbl,euDis,dataForDisplay.rangeNames{i},refType,useMedianFlag);
end

%%  To find powermatched mid and old list
% powerMatchedSubjectNameLists = getPowerMatchedSubjectList(folderSourceString,subjectNameListFinal,projectName,refType,protocolType);