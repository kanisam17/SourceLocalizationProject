% prepareLORETAFiles
% Requires a 1x2 cell array of subjectNameListFinal
% input from the result of "ConvertMATtoDAT.m" code
clc;clear;
load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\subjectNameListFinal.mat');
dataStr{1} = 'mid'; dataStr{2} = 'old';
% subjectNameListFinal = load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\ADGammaProjectCodes\subjectNameListFinal.mat');
%subjectNameListFinal = caseList.subjectNameListMatched;
folderSourceString = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10';
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP

% LORETA files
folderLORETA = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10\interpolatedData\saved_data\text';
folderName = 'newData';


x = load('goodProtFlag');

for i=1:size(dataStr,2)
    subjectNames = subjectNameListFinal{i};  %default is subjectNameListFinal{i};
    
    folderToSaveBL = fullfile(folderLORETA,folderName,dataStr{i},'BL');
    folderToSaveST = fullfile(folderLORETA,folderName,dataStr{i},'ST');
    mkdir(folderToSaveBL); mkdir(folderToSaveST);
    
    %     for j=1:length(subjectNames)
    %         subjectName = subjectNames{j};
    %         subjectName = subjectNames{j};%added by kan
    
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        %         for iSub = 1:length(subjectNames{1,j})%added by kan
        %             subjectNametemp = subjectNames{j};
        %             subjectName = subjectNametemp{iSub};%added by kan
        
        [expDates,protocolNames] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);
        
        pos = find(strcmp(x.uniqueSubjectNames,subjectName));
        goodProts = x.goodProtFlagList{pos};
        
        for k=1:length(protocolNames)
            if goodProts(k) % If this is a good protocol, % Retrieve filename from decimated data and save in the output folder
                data = load (fullfile(folderSourceString,projectName,protocolType,...
                    [subjectName '-' expDates{k} '-' protocolNames{k} '.mat']));
                %saveBL
                bl = data.eegData(:,:,126:250); %for BL find(timeVals>-0.5 & timeVals<0)
                a = permute(bl,[1,3,2]);
                y = reshape (a,[size(a,1),size(a,2)*size(a,3)]);
                writematrix(y', fullfile(folderToSaveBL,[subjectName '-' expDates{k} '-' protocolNames{k} '.txt']), 'Delimiter' , 'space');
                %saveST
                st = data.eegData(:,:,314:438); % for ST find(timeVals>0.25 & timeVals<0.75)
                a1 = permute(st,[1,3,2]);
                y1 = reshape (a1,[size(a,1),size(a,2)*size(a,3)]);
                writematrix(y1', fullfile(folderToSaveST,[subjectName '-' expDates{k} '-' protocolNames{k} '.txt']), 'Delimiter' , 'space');
            end
        end
    end
    %end
end

