function matchedSubjectNameLists = getPowerMatchedSubjectList(folderSourceString,subjectNameListFinal,projectName,refType,protocolType,freqBand)

powerMatchingBinWidth = 0.5; % dB
dataForDisplay{1} = combineAnalyzedData(folderSourceString,subjectNameListFinal{1},projectName,refType,protocolType);
dataForDisplay{2} = combineAnalyzedData(folderSourceString,subjectNameListFinal{2},projectName,refType,protocolType);
freqBand = 3; %change frequency
tmp = (dataForDisplay{1}.powerDBAllSubjects); dataToMatch{1} = tmp(:,freqBand);
tmp = (dataForDisplay{2}.powerDBAllSubjects); dataToMatch{2} = tmp(:,freqBand);
minVal = min(min(dataToMatch{1}),min(dataToMatch{2}));
maxVal = max(max(dataToMatch{1}),max(dataToMatch{2}));
powerBins = minVal:powerMatchingBinWidth:maxVal;
numPowerBins = length(powerBins);
d = powerBins(2)-powerBins(1);
powerBins = [powerBins powerBins(numPowerBins)+d];

matchedSubjectNameLists{1} = [];
matchedSubjectNameLists{2} = [];

for i=1:numPowerBins
    pos1 = intersect(find(dataToMatch{1}>=powerBins(i)),find(dataToMatch{1}<powerBins(i+1)));
    pos2 = intersect(find(dataToMatch{2}>=powerBins(i)),find(dataToMatch{2}<powerBins(i+1)));
    [equalPos1,equalPos2] = getEqualNumOfIndices(pos1,pos2);
        matchedSubjectNameLists{1} = cat(2,matchedSubjectNameLists{1},subjectNameListFinal{1}(equalPos1));
        matchedSubjectNameLists{2} = cat(2,matchedSubjectNameLists{2},subjectNameListFinal{2}(equalPos2));
end  
end

function [x2,y2] = getEqualNumOfIndices(x1,y1)

N1 = length(x1);
N2 = length(y1);

if (N1==0) || (N2==0) % one of the two is an empty array
    x2=[]; y2=[];
    
elseif N1==N2
    x2=x1; y2=y1;
    
elseif N1<N2
    x2 = x1;
    randVals = randperm(N2);
    y2 = y1(sort(randVals(1:N1)));
    
else %N1>N2
    y2 = y1;
    randVals = randperm(N1);
    x2 = x1(sort(randVals(1:N2)));
end
end
