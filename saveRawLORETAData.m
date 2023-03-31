% saveRawLORETAData
% Retrives single trial sLORETA data from the relevant folders and saves
% the following:
% 1. mean (across trials) power for BL and ST
% t-stats and p-vals after performing t-test.

function saveRawLORETAData(subjectNameListFinal,dataStr,folderLORETA,folderOutput)

% x = load('goodProtFlag');
x = load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\goodProtFlag');

for i=1:2
    subjectNames = subjectNameListFinal{i};

    folderToSave = fullfile(folderOutput,dataStr{i});
    makeDirectory(folderToSave);
    
    folderLORETABL = fullfile(folderLORETA,dataStr{i},'BL');
    folderLORETAST = fullfile(folderLORETA,dataStr{i},'ST');
    
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
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
%                 %Deltapower
%                 
%                 DeltaP = 10*(log10(sourceData(i).ST) - log10(sourceData(i).BL));
%                 
            end
        end

        
        clear mDataBL mDataST tStats pVals numStimuli
        numStimuli = size(dataBLtmp,1);
        mDataBL = squeeze(mean(dataBLtmp,1));
        mDataST = squeeze(mean(dataSTtmp,1));
%         mDataDeltaP = squeeze(mean(dataDeltaP,1));
        
        % do statistical testing
        numFreqRanges = size(dataBLtmp,2);
        numVoxels = size(dataBLtmp,3);
        tStats = zeros(numFreqRanges,numVoxels);
        pVals = zeros(numFreqRanges,numVoxels);
        dataDeltaP = zeros(numFreqRanges,numVoxels);
        for nF=1:numFreqRanges
            for nV=1:numVoxels
                [~,pVals(nF,nV),~,stats] = ttest(squeeze(dataSTtmp(:,nF,nV)),squeeze(dataBLtmp(:,nF,nV)));
                tStats(nF,nV) = stats.tstat;
                %save deltaPower
%                 dataDeltaP(nF,nV)= (squeeze(10*(log10((dataSTtmp(:,nF,nV))) - log10((dataBLtmp(:,nF,nV))))))';
                
            end
        end
        
        % Save data
        fileToSave = fullfile(folderToSave,subjectName);
        save(fileToSave,'mDataBL','mDataST','tStats','pVals','numStimuli');
    end
end