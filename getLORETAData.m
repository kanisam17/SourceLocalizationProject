% getLORETAData
% Retrives sLORETA data from the relevant folder.

function [allDataBL,allDataST] = getLORETAData(subjectNameListFinal,dataStr,folderLORETA)

x = load('goodProtFlag');

allDataBL = cell(1,2);
allDataST = cell(1,2);

for i=1:2
    disp(dataStr{i});
    subjectNames = subjectNameListFinal{i};

    folderToSaveBL = fullfile(folderLORETA,dataStr{i},'BL','Averaged','sLORtoText');
    folderToSaveST = fullfile(folderLORETA,dataStr{i},'ST','Averaged','sLORtoText');
    
    clear dataBL dataST
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        
        [expDates,protocolNames] = getProtocolDetailsForAnalysis('ADGammaProject',subjectName,'SF_ORI');
        
        pos = strcmp(x.uniqueSubjectNames,subjectName);
        goodProts = x.goodProtFlagList{pos};
        
        count=1; 
        clear dataBLtmp dataSTtmp
        for k=1:length(protocolNames)   
            if goodProts(k) % If this is a good protocol
                fileNameBL = fullfile(folderToSaveBL,[subjectName '-' expDates{k} '-' protocolNames{k} '-slor.txt']);
                dataBLtmp(count,:,:) = textread(fileNameBL); %#ok<*AGROW,*DTXTRD>
                fileNameST = fullfile(folderToSaveST,[subjectName '-' expDates{k} '-' protocolNames{k} '-slor.txt']);
                dataSTtmp(count,:,:) = textread(fileNameST);
                count=count+1;
            end
        end
        
        dataBL(j,:,:) = squeeze(mean(dataBLtmp,1));
        dataST(j,:,:) = squeeze(mean(dataSTtmp,1));
    end
    allDataBL{i} = dataBL;
    allDataST{i} = dataST;
end
