% saveRawLORETAData
% Retrives single trial sLORETA data from the relevant folders and saves
% the following:
% 1. mean (across trials) power for BL and ST
% t-stats and p-vals after performing t-test.

function saveRawLORETArandTrialData(subjectNameListFinal,dataStr,folderLORETA,folderOutput,trialIdxListFinal)

% x = load('goodProtFlag');
x = load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\goodProtFlag');

for i=1:length(dataStr)
    
    subjectNames = subjectNameListFinal{i};
    
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        
        folderToSave = fullfile(folderOutput,dataStr{i});
        
        folderLORETABL = fullfile(folderLORETA,dataStr{i},'BL');
        folderLORETAST = fullfile(folderLORETA,dataStr{i},'ST');
        
        %     for j=1:length(subjectNames) % for case and control averaged
        %         subjectName = subjectNames{j};
        
        
        
        disp([ dataStr{i} ', ' subjectName]);
        
        [expDates,protocolNames] = getProtocolDetailsForAnalysis('ADGammaProject',subjectName,'SF_ORI');
        
        pos = strcmp(x.uniqueSubjectNames,subjectName);
        goodProts = x.goodProtFlagList{pos};
        
        dataBLtmp = [];
        dataSTtmp = [];
        dataDeltaP = []; %Added by Kan
        
        for k=1:length(protocolNames)
            if goodProts(k) % If this is a good protocol
                %                 filesBL = dir(fullfile(folderLORETABL,[subjectName '-' expDates{k} '-' protocolNames{k}],'sLORtoText\*.txt'));
                folder2ReadBL = fullfile(folderLORETABL,[subjectName '-' expDates{k} '-' protocolNames{k}],'sLORtoText');
                for l = 1:size(trialIdxListFinal{i}{j}{k},2)
                    epochNum = trialIdxListFinal{i}{j}{k}(l);
                    file2ReadBL = [subjectName '-' expDates{k} '-' protocolNames{k} '-epoch' num2str(epochNum,'%04d') '-slor.txt'];
                    tmpBL(l,:,:) = textread(fullfile(folder2ReadBL,file2ReadBL));
                end
            
                % clear tmpBL
                %                 for f=1:size(filesBL,1)
                %                     tmpBL(f,:,:) = textread(fullfile(filesBL(f).folder,filesBL(f).name)); %#ok<*DTXTRD,*AGROW>
                %                 end
                dataBLtmp = cat(1,dataBLtmp,tmpBL);
                
                
                folder2ReadST = fullfile(folderLORETAST,[subjectName '-' expDates{k} '-' protocolNames{k}],'sLORtoText');
                for l = 1:size(trialIdxListFinal{i}{j}{k},2)
                    epochNum = trialIdxListFinal{i}{j}{k}(l);
                    file2ReadST = [subjectName '-' expDates{k} '-' protocolNames{k} '-epoch' num2str(epochNum,'%04d') '-slor.txt'];
                    tmpST(l,:,:) = textread(fullfile(folder2ReadST,file2ReadST));
                end
                
%                 filesST = dir(fullfile(folderLORETAST,[subjectName '-' expDates{k} '-' protocolNames{k}],'sLORtoText\*.txt'));
%                 
%                 clear tmpST
%                 for f=1:size(filesST,1)
%                     tmpST(f,:,:) = textread(fullfile(filesST(f).folder,filesST(f).name));
%                 end
                dataSTtmp = cat(1,dataSTtmp,tmpST);
                
                dataDeltaP = 10*(log10(dataSTtmp) - log10(dataBLtmp)); %Added by Kan 20/02/2023
            end
        end
        
        clear mDataBL mDataST tStats pVals numStimuli
        numStimuli = size(dataBLtmp,1);
        mDataBL = squeeze(mean(dataBLtmp,1));
        mDataST = squeeze(mean(dataSTtmp,1));
        mDataDeltaP = squeeze(mean(dataDeltaPtmp,1)); % Added by kan
        
        
        % do statistical testing
        numFreqRanges = size(dataBLtmp,2);
        numVoxels = size(dataBLtmp,3);
        tStats = zeros(numFreqRanges,numVoxels);
        pVals = zeros(numFreqRanges,numVoxels);
        for nF=1:numFreqRanges
            for nV=1:numVoxels
                [~,pVals(nF,nV),~,stats] = ttest(squeeze(dataSTtmp(:,nF,nV)),squeeze(dataBLtmp(:,nF,nV)));
                tStats(nF,nV) = stats.tstat;
            end
        end
        
        % Save data
        mkdir (folderToSave); 
        fileToSave = fullfile(folderToSave,subjectName); 
        save(fileToSave,'mDataBL','mDataST','tStats','pVals','numStimuli');
    end
end
end

%% Average of matched subject w.r.t. case subjects %added by kan
% folderSource = fullfile('N:\Projects\Kanishka_SourceLocalizationProject\data\sLORETA_Thres10\caseControl\control');
%
% load('caseListAgeMatched.mat')
%
%
%
% for i = 1:size(caseList.subjectNameListMatched{1},2)
%     for j = 1:size(caseList.subjectNameListMatched{1, 1}{1,i},2)
%         for k = 1:size(caseList.subjectNameListMatched{1, 1}{1,i}{1,j})
%         data = load(caseList.subjectNameListMatched{1, 1}{1,i}{1,j});
%         mDataBL_mean = zeros(3,6239,size(caseList.subjectNameListMatched{1, 1}{1,i},2));
%         mDataST_mean = zeros(3,6239,size(caseList.subjectNameListMatched{1, 1}{1,i},2));
%         pVals_mean = zeros(3,6239,size(caseList.subjectNameListMatched{1, 1}{1,i},2));
%         tStats_mean = zeros(3,6239,size(caseList.subjectNameListMatched{1, 1}{1,i},2));
%
%
%         end
%     end
% end




