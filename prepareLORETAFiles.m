% sLORETA is run using a software: {Kanishka to provide more details}
% This program takes decimatedData and writes it in a format that is used
% by the software. The output of the software is in the same folder where
% input files are kept.

% Requires a 1x2 cell array of subjectNameListFinal and another array
% dataStr that provides the labels of the groups.

function prepareLORETAFiles(subjectNameListFinal,dataStr,folderLORETA)

if ~exist('dataStr','var');         dataStr=[];                         end

if isempty(dataStr)
    dataStr{1} = 'mid'; 
    dataStr{2} = 'old';
end

folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject\decimatedData';
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP

x = load('goodProtFlag');

for i=1:2
    subjectNames = subjectNameListFinal{i};

    folderToSaveBL = fullfile(folderLORETA,dataStr{i},'BL');
    folderToSaveST = fullfile(folderLORETA,dataStr{i},'ST');
    
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        
        [expDates,protocolNames] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);
        
        pos = find(strcmp(x.uniqueSubjectNames,subjectName));
        goodProts = x.goodProtFlagList{pos};
        
        for k=1:length(protocolNames)
            if goodProts(k) % If this is a good protocol
                % Retrieve filename from decimated data and save in the
                % output folder - Kanishka to update
            end
        end
    end
end
end
