% getLORETAData
% Retrives sLORETA data from the relevant folder.

function [allDataBL,allDataST,allDataDeltaP,alltStats,allpVals] = getLORETAData(subjectNameListFinal,dataStr,folderLORETA)
numGroups = length(subjectNameListFinal);
allDataBL = cell(1,numGroups);
allDataST = cell(1,numGroups);
allDataDeltaP = cell(1,numGroups);
alltStats = cell(1,numGroups);
allpVals = cell(1,numGroups);

for i=1:numGroups
    subjectNames = subjectNameListFinal{i};
    clear dataBL dataST tStats pVals
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        tmp = load(fullfile(folderLORETA,dataStr{i},[subjectName '.mat']));
        dataBL(:,:,j) = tmp.mDataBL; %#ok<*AGROW>
        dataST(:,:,j) = tmp.mDataST;
        dataDeltaP(:,:,j) = 10*(log10(tmp.mDataST) - log10(tmp.mDataBL));
        tStats(:,:,j) = tmp.tStats;
        pVals(:,:,j) = tmp.pVals;
    end
    allDataBL{i} = dataBL;
    allDataST{i} = dataST;
    allDataDeltaP{i} = dataDeltaP;
    alltStats{i} = tStats;
    allpVals{i} = pVals;
end

