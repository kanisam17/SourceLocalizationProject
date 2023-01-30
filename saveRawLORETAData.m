% saveRawLORETAData
% Retrives single trial sLORETA data from the relevant folders and saves
% the following:
% 1. mean (across trials) power for BL and ST
% t-stats and p-vals after performing t-test.

function saveRawLORETAData(subjectNameListFinal,dataStr,folderLORETA,folderOutput)

x = load('goodProtFlag');

for i=1:length(dataStr)  
  
    subjectNames = subjectNameListFinal{i}; 

    for j=1:length(subjectNames)
        for iSub = 1:length(subjectNames{1,j})%added by kan
            subjectNametemp = subjectNames{j};
            subjectName = subjectNametemp{iSub};%added by kan
    
    folderToSave = fullfile(folderOutput,dataStr{i});
    makeDirectory(folderToSave);
    
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
        
        for k=1:length(protocolNames)   
            if goodProts(k) % If this is a good protocol
                % Read all files in these folders
                filesBL = dir(fullfile(folderLORETABL,[subjectName '-' expDates{k} '-' protocolNames{k}],'sLORtoText\*.txt'));
                clear tmpBL
                for f=1:size(filesBL,1)  
                    tmpBL(f,:,:) = textread(fullfile(filesBL(f).folder,filesBL(f).name)); %#ok<*DTXTRD,*AGROW>
                end
                dataBLtmp = cat(1,dataBLtmp,tmpBL);
                
                filesST = dir(fullfile(folderLORETAST,[subjectName '-' expDates{k} '-' protocolNames{k}],'sLORtoText\*.txt'));
                clear tmpST
                for f=1:size(filesST,1)  
                    tmpST(f,:,:) = textread(fullfile(filesST(f).folder,filesST(f).name));
                end
                dataSTtmp = cat(1,dataSTtmp,tmpST);
            end
        end

        clear mDataBL mDataST tStats pVals numStimuli
        numStimuli = size(dataBLtmp,1);
        mDataBL = squeeze(mean(dataBLtmp,1));
        mDataST = squeeze(mean(dataSTtmp,1));
       
        
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
        
        
        
        
        