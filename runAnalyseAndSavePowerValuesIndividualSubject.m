% Modified from runAnalyseAndSaveValuesIndividualSubject in ageProject.
% Here, we simply save the power values, which are subsequently needed to
% select subjects for sourceLocalization analysis.

% Data avalable in the decimatedData folder needs to be properly extracted
% for analysis. This program extracts the part that is used for final
% analysis and saves it under analyzedData.

clear; clc;

% Mandatory fixed options
folderSourceString = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10\interpolatedData'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];


% Choose one of these options
refType = 'unipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjects = getGoodSubjectsProjectwise(projectName,1);
uniqueSubjectNames = getGoodFileNamesForSubjects(goodSubjects{1});

for iSub = 1:length(uniqueSubjectNames)
    subjectName = uniqueSubjectNames{iSub};
    disp([num2str(iSub) ': ' subjectName]);
    analyseAndSaveValuesIndividualSubject(folderSourceString,subjectName,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag); % Save data in analyzedData
end
