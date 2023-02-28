%figure for euclidian distance binned; % remove if nesting in other code #kan
% % % % % % plotBinWidth = 35; plotBinLimit = 141;
% % % % % % % freqRange = 3; % choose the frequency range
% % % % % % useMedianFlagData = true;
% % % % % % vrbl = dataDeltaP; % dataDeltaP tStats change variable
% % % % % % for i = 1:size(strList,2) %group
% % % % % %     for j = 1:size(vrbl{i},3) % sub
% % % % % %         for k = 1:numFreqRanges% freq
% % % % % %         mtStats= vrbl{i}(k,:,j);
% % % % % %         meuDis = euDis{i}(k,:,j);
% % % % % %        [~,mean_binY1{i}(k,:,j)] = plotIndivConnData(mtStats,meuDis,plotBinWidth,plotBinLimit,[0 0 0],useMedianFlagData,0);
% % % % % %         end
% % % % % %     end
% % % % % % end

%%
% Preallocate output cell array
mtStats = cell(1, 2);
meuDis = cell(1,2);
% Loop over each cell of the input tStats array
vrbl = dataDeltaP; % dataDeltaP or tStats change variable
for i = 1:numel(tStats)
    % Calculate the mean over the third dimension (subjects)
    meanTStats = mean(vrbl{i}, 3);
    meanEuDis = mean (euDis{i},3); 
    % Reshape the result to be numFreq x numVoxels
    mtStats{i} = reshape(meanTStats, size(meanTStats, 1), []);
    meuDis{i} = reshape(meanEuDis, size(meanEuDis, 1), []);
    
    
end
%% plot figures 
dataForStats = cell(1,3);
plotBinWidth = 14; plotBinLimit = 141;


for numFreq = 1:3
    figure;
    for numGroup = 1:2
        colorMrkr = {[0 0 0],[0.5 0.5 0.5]};
        alphaErrorbar = 0.65; % transparency level for errorbars
        headingColor = colorMrkr;
        useMedianFlagData = true;
        [h{numGroup},binned_Y{numGroup}] = plotIndivConnData(tStats{numGroup}(numFreq,:,:),euDis{numGroup}(numFreq,:,:),plotBinWidth,plotBinLimit,colorMrkr{numGroup},useMedianFlagData,1);
        dataForStats{numFreq} = binned_Y;
        hold on
    end
    
end

%% statistics
% Set the confidence level
alpha = 0.05;
% Get the number of frequency bands and bins
[numFreq, numBins, ~] = size(dataForStats{1});
% Initialize the p-value matrix
p_mat = NaN(3, 10);%p_mat = NaN(numFreq, numBins);
% Loop through each frequency band and bin
for i = 1:3%numFreq
    for j = 1:10%numBins
        data1 = squeeze(dataForStats{i}{1,1}(:,j));
        data2 = squeeze(dataForStats{i}{1,2}(:,j));
        % Concatenate the data 
        data = [data1; data2];
        group = [ones(size(data1)); 2*ones(size(data2))];
        % Perform the Kruskal-Wallis test
        [p, ~, ~] = kruskalwallis(data, group, 'off');
        % Store the p-value in the matrix
        p_mat(i,j) = p;
    end
end
%% 

function [h,mean_binY1] = plotIndivConnData(mtStats,meuDis,binWidth,binLimit,colorMrkr,medianFlag,plotSwitch)
%% variables list
% tStats = subjectwise t values of each voxel, subjects: 76, voxels:6239, numFreq=3, size: (3 x 6239 x 76)
% euDis = subjectwise euc. distance from max activ voxel to each voxel, dimension same as above
% binWidth = 10 as max value for euDis is 150, binned in 10 bins
% colorMrkr = used for two different groups
% plotSwitch = label to either plot or not
if ~exist('plotSwitch','var');    plotSwitch = 1;       end
connectDiscrete =  '-o';
DiscreteVisibility = 'on';
binEdges = 0:binWidth:binLimit;
nbins = length(binEdges)-1;
binned_bin_meuDis = binEdges(1:end-1)+(binWidth/2);

nsub = size(mtStats,3);
% nsub = 3; % 
% nvoxels = size(tStats,1);

mean_binY1 = zeros(nsub,nbins); %%default mean_binY1 = zeros(nsubjects,nbins);

% median_binned_binY_plot = zeros(1,nbins);%%%%% save this variable
% std_binned_binY_plot = zeros(1,nbins);%%%%% save this variable




 for i = 1:nsub
    binned_meuDis = discretize(meuDis(:,:,i),binEdges); %%binned_meuDis = discretize(meuDis(:,i),binEdges);
    binned_mtStats = cell(1,nbins); 
    test_mtStats = mtStats(:,:,i);

    for b = 1:nbins % number of bins
        binned_mtStats{b} = test_mtStats(binned_meuDis == b);
    end


    if(medianFlag)
        %     median_binned_binY = nanmedian(mean_binY1,1);
        mSEM = @(data)std(bootstrp(1000,@nanmedian,data));
    else
        %     median_binned_binY = nanmean(mean_binY1,1);
        mSEM = @(data)std(bootstrp(1000,@nanmean,data));
    end
    if(medianFlag)
        mean_binY1(i,:) = cellfun(@nanmedian,binned_mtStats); % Note:add if bin value is empty,dont assign SEM.
        std_binned_binY(i,:)= cellfun(mSEM,binned_mtStats);
        
    else
        mean_binY1(i,:) = cellfun(@nanmean,binned_mtStats);
        std_binned_binY(i,:)= cellfun(mSEM,binned_mtStats);
       
    end
    
    median_binned_binY(i,:) = mean_binY1(i,:);
    %
    % std_binned_binY = mSEM(mean_binY1);
 end
  disp(size(std_binned_binY))
  std_binned_binY_plot = mean(std_binned_binY);
  median_binned_binY_plot = mean(median_binned_binY);
 
 if(plotSwitch)
        if(length(std_binned_binY_plot)==1)
            h = errorbar(binned_bin_meuDis,median_binned_binY_plot,nan(1,8),connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
            legend('mid','old'); xlabel('eucDistance(mm)'); ylabel('deltaP')
        else
            h = errorbar(binned_bin_meuDis,median_binned_binY_plot,std_binned_binY_plot,connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
            legend('mid','old'); xlabel('eucDistance(mm)'); ylabel('deltaP')
        end
    else
        h = [];
    end
end