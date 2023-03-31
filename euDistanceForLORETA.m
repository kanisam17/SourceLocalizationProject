function [euDis,tStats,dataDeltaP]= euDistanceForLORETA(folderLORETA,subjectNameListFinal,strList,xyz,useMedianFlag,idxFPTO)
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
% clear;clc
% folderLORETA = 'N:\Projects\Kanishka_SourceLocalizationProject\data\sLORETA_Thres10';
% numGroups = {'mid','old'};
% xyz = xlsread ('voxelInfo.xlsx');%}

%% Subject Listing 1. All Subject thresh-10
% for i = 1:size(strList,2)
%     files{i} = dir(fullfile(folderLORETA,strList{i},'*.mat'));
%     for j = 1:size(files{i},1)
%         temp = load (fullfile(files{i}(j).folder, files{i}(j).name), 'tStats');
%         temp1(:,:,j) = temp.tStats;
%         end
%         tStats{i} = temp1;
%
% end
clear tStats; clear euDis;
% %for Subject list matched
for i = 1:length(strList)
    temp{i}=  subjectNameListFinal{i};
    for j = 1:length(temp{i})
        g{i}{j} = [temp{i}{j} '.mat'];
        tempM = load(fullfile(folderLORETA,strList{i},g{i}{j}),'tStats');
        temp1(:,:,j) = tempM.tStats;
        
        tempBL = load(fullfile(folderLORETA,strList{i},g{i}{j}),'mDataBL');%added by Kan 20/02/2023
        tempST = load(fullfile(folderLORETA,strList{i},g{i}{j}),'mDataST');%added by Kan 20/02/2023
        deltaP(:,:,j) = 10*(log10(tempST.mDataST) - log10(tempBL.mDataBL));%added by Kan 20/02/2023
    end
    tStats{i} = temp1;
    dataDeltaP{i} = deltaP;
    %
    %     if useMedianFlag
    %         tStats{i} = median(temp1,3);
    %     else
    %         tStats{i} = mean(temp1,3);
    %      end
    
end


%% calualate euDis
vrbl = dataDeltaP; % dataDeltaP/tStats, change the vrbl to define euDis.
for i = 1:length(strList)
    if useMedianFlag
        tempV{i} = median(vrbl{i}, 3);
    else
        tempV{i} = mean(vrbl{i}, 3);
    end
end

for i = 1:length(strList)
    tempV{i}(idxFPTO) = NaN; %% make all sub lobar, limbic,frontal(optional) values as Nan
%     tempV{i}(idxT) = NaN; %% make all temporal(optional) values as Nan
%     tempV{i}(idxF) = NaN; %% make all frontal(optional) values as Nan
end

%% Concatenating method for euDis
catTempV = cat(2,tempV{:}); % concatenation of all subjects data
[~, maxPosition] = (max(catTempV,[],2));
%position=squeeze(position);
[~, minPosition] = (min(catTempV,[],2));
position = [minPosition(1);maxPosition(2:3)];
reptxyz = [xyz;xyz];

%%

for i = 1:length(strList)
    for j = 1:3%(freqType)
        for k = 1:length(subjectNameListFinal{i}) %%%%%%temp
            a = position(j,:);
            
            rept = repmat(reptxyz(a,:),[6239,1]); % create matrix of repeted xyz coordinates
            euDis{i}(j,:,k)= sqrt(sum((rept - xyz).^ 2,2));
        end
    end
end




% plot euDis to tStats
%average for groups
z = mean(tStats{1,1},3);
z1 = mean(tStats{1,2},3);
u = mean(euDis{1,1},3);
u1 = mean(euDis{1,2},3);

% mtStats{1} = z; mtStats{2} = z1;
% meuDis{1} =  u; meuDis{2} = u1;
clear z z1 u u1;

% % % % % %% Plot%%
% % for i=1:2%numGroups
% %      for j = 1%:3%numFreqRanges
% %          scatter(meuDis{1}(j,:), mtStats{1}(j,:),'r','Parent',hPlots(j,9));hold on
% %          scatter(meuDis{2}(j,:), mtStats{2}(j,:),'k','Parent',hPlots(j,9));
% %      end
% %  end

