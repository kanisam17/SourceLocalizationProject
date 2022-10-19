% getLORETAData
% Retrives sLORETA data from the relevant folder.

function [allDataBL,allDataST,alltStats,allpVals] = getLORETAData(subjectNameListFinal,dataStr,folderLORETA)

allDataBL = cell(1,2);
allDataST = cell(1,2);
alltStats = cell(1,2);
allpVals = cell(1,2);

for i=1:2
    subjectNames = subjectNameListFinal{i};

    clear dataBL dataST tStats pVals
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        tmp = load(fullfile(folderLORETA,dataStr{i},[subjectName '.mat']));
        
        dataBL(j,:,:) = tmp.mDataBL; %#ok<*AGROW>
        dataST(j,:,:) = tmp.mDataST;
        tStats(j,:,:) = tmp.tStats;
        pVals(j,:,:) = tmp.pVals;
    end
    allDataBL{i} = dataBL;
    allDataST{i} = dataST;
    alltStats{i} = tStats;
    allpVals{i} = pVals;
end
