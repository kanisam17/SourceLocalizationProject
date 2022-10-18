% Returns the indices for areas listed in areaList. The original
% list has 6 areas - Occipital, Parietal, Temporal, Frontal, Limbic and
% Sub-lobar

function [posList,xyz,areaList] = getVoxelInfo(areaList,displayFlag)

if ~exist('areaList','var');            areaList=[];                    end
if ~exist('displayFlag','var');         displayFlag=0;                  end

if isempty(areaList)
    areaList{1} = 'Occipital';
    areaList{2} = 'Parietal';
    areaList{3} = 'Temporal';
    areaList{4} = 'Frontal';
end

[allCoordinates,tmp] = xlsread('voxelInfo.xlsx');
xyz = allCoordinates;
labels = tmp(:,1);
uniqueLabels = unique(labels);
numUniqueLabels = length(uniqueLabels);

allLabelNums = zeros(1,size(xyz,1));
for i=1:length(uniqueLabels)
     allLabelNums(strcmp(labels,uniqueLabels{i}))=i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numAreas = length(areaList);
posList = cell(1,numAreas);
colorList = jet(numAreas);

for i=1:numAreas
    tmp=strfind(uniqueLabels,areaList{i});
    matchArray = zeros(1,numUniqueLabels);
    for j=1:numUniqueLabels
        matchArray(j) = ~isempty(tmp{j});
    end
    goodLabelNum = find(matchArray);
    posList{i} = find(allLabelNums==goodLabelNum);
    
    if displayFlag
        p=posList{i}; 
        scatter3(xyz(p,1),xyz(p,2),xyz(p,3),100,colorList(i,:));
        hold on;
    end
end
end