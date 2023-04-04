%% Find total number of stimuli/trials available for each subject. Choose ones 
function listNumStimuli = trialsAvailableEachSubject()
% clear;clc;
%folderLORETA = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\sourceData\LORETA\data\Age';
folderLORETA = 'N:\Projects\Kanishka_SourceLocalizationProject\data\interpolatedData';
numGroup = {'mid','old'};

for i = 1:size(numGroup,2)
    files{i} = dir(fullfile(folderLORETA,numGroup{i},'*.mat'));
    for j = 1:size(files{i},1)
        temp = load (fullfile(files{i}(j).folder, files{i}(j).name), 'numStimuli');
        temp1(j) = temp.numStimuli;
    end
    listNumStimuli{i} = temp1; % Trials available for each subject
    maxStim{i} = max(listNumStimuli{i});
    minStim{i} = min(listNumStimuli{i});
    
end
%figure; plot (listNumStimuli{1, 1});hold on; plot (listNumStimuli{1, 2});
% figure;  
% plot (listNumStimuli{1, 1}(listNumStimuli{1, 1}>180));
% hold on;
% plot (listNumStimuli{1, 2}(listNumStimuli{1, 2}>180));
% 
% 
% sum(listNumStimuli{1, 1}>200);
% sum(listNumStimuli{1, 2}>200);
% 
% trialMatch1Pos = (listNumStimuli{1}>180) & ageGroup1Pos; strList{1} = 'mid';
% trialMatch2Pos = (listNumStimuli{2}>180) & ageGroup2Pos; strList{2} = 'old';
% 
% powerMatchedSubjectPos = ageGroup1Pos | ageGroup2Pos; % use when ageGroup1Pos and ageGroup1Pos are generated using power matched List
% 
% 
% y1 = ismember(uniqueSubjectNames, subjectNameListFinal{1});
% y2 = ismember(uniqueSubjectNames, subjectNameListFinal{2});
% y = y1 | y2; % y is the logical index of all mid and old sub togather

