% prepareLORETAFiles
% Requires a 1x2 cell array of subjectNameListFinal
dataStr{1} = 'mid'; dataStr{2} = 'old';

folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject\decimatedData';
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP

% LORETA files
folderLORETA = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\Kanishka_SourceLocalizationProject\data';
folderName = 'sLORETA_Thres0';

x = load('goodProtFlag');

for i=1:2
    subjectNames = subjectNameListFinal{i};

    folderToSaveBL = fullfile(folderLORETA,folderName,dataStr{i},'BL');
    folderToSaveST = fullfile(folderLORETA,folderName,dataStr{i},'ST');
    
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        
        [expDates,protocolNames] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);
        
        pos = find(strcmp(x.uniqueSubjectNames,subjectName));
        goodProts = x.goodProtFlagList{pos};
        
        for k=1:length(protocolNames)
            if goodProts(k) % If this is a good protocol
                % Retrieve filename from decimated data and save in the output folder
            end
        end
    end
end
