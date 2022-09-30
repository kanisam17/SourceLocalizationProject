% For each subject, bad protocols are rejected but this information was
% initially not saved anywehre. But later it is saved in the connectivityProject. This
% program simply reads this information from there and saves it locally.

% Needs to be run only once.

folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP

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

folderName = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject\TLSAEEGProjectPrograms\connectivityProjectCodes\analyzedData\ConnectivityProject\SF_ORI';
for i=1:length(uniqueSubjectNames)
    x = load(fullfile(folderName,[uniqueSubjectNames{i} '_unipolar_stRange_250_750_RemoveSF1_ppc.mat']));
    goodProtFlagList{i} = x.goodProtFlag; %#ok<*SAGROW>
    numGoodTrialsList{i} = x.numGoodTrials;
end

save('goodProtFlag.mat','uniqueSubjectNames','goodProtFlagList','numGoodTrialsList');