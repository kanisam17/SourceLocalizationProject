% temp code to get trialIdx of all subjects for all protocols 
% run this code only after running the runAnalyseAndSaveValuesIndividualSubject.m


function [trialIdxAllSub,subjectNamesOnly,subNameIdx,goodProtFlag] = getTrialIdxAndGoodSubject64Elecs(folderSourceString,projectName,protocolType,uniqueSubjectNames)

files = dir(fullfile(folderSourceString, projectName, protocolType, '*.mat')); 
% files= files(3:end);% to select only subject mat files

% 

for i = 1:length(files)
    file_lengths(i) = length(files(i).name);
end
files = files(file_lengths<42); % 42 used to get only files with filename as '005_SR_F2_unipolar_stRange_250_750.mat'. Where 42 is the max 
...character length to find this type of filename  
    
filenames = {files.name}; % extract file names
index = strfind(filenames, 'unipolar'); % find index of 'unipolar' in filenames

%%

for i = 1:length(files)
    data = load(fullfile(files(i).folder, files(i).name));
    %     for j = 1:length(data.trialIdx)
    try
        trialIdxAllSub{i} = data.trialIdx;
    catch
    end
    subjectNamesOnly{i} = filenames{i}(1:(index{i}-2));
    goodProtFlag{i} = data.allProtocolsBLData.goodProtFlag;  
end
subNameIdx = ismember(subjectNamesOnly,uniqueSubjectNames);





