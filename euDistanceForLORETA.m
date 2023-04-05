function euDis = euDistanceForLORETA(vrbl,xyz,useMedianFlag)
%{ 1. Find euclidian distance 'euDis' from the highest activated voxel to all
% other voxels.
% 2. The euDis = sqrt(sum((A - B) .^ 2)) where A is the max activated voxel
% coordinates(X0,Y0,Z0) and B is corrdinates of other voxel(X1,Y1,Z1).
% 3. Plot scatter for 'euDis' on x axis and CSD/Power value for voxels on Y axis.
% Prerequisite: matrix of (3 x 6239) for each group (Mid/Old, Control/Case,
% MatchedPower Mid/Old)for CSD/Power value and (3 x 6239) for euDis
% 4. Make average of points on x axis to reduce cluster.
% For euDis, first take the max of tStats and then find position of that
% voxel, then calcualte the euDis.

numGrps = length(vrbl);
numFreq = size(vrbl{1},1);
numVoxel = size(vrbl{1},2);
euDis = zeros(numFreq,numVoxel);

% calualate euDis
accumVrbl = cat(3,vrbl{1},vrbl{2});

if useMedianFlag
    tempVavg = median(accumVrbl,3);
else
    tempVavg = mean(accumVrbl,3);
end

[~, minPosition] = (min(tempVavg(1,:),[],2));
[~, maxPosition] = (max(tempVavg(2:3,:),[],2));
position = [minPosition(1);maxPosition(1:2)];

for j = 1:numFreq
%     rept = repmat(xyz(position(j),:),[numVoxel,1]); % create matrix of repeted xyz coordinates
    euDis(j,:)= sqrt(sum((xyz - xyz(position(j),:)).^ 2,2));
end
end

