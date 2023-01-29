
clear;clc;
%%
% Run in data (decimated) folder
files = dir('*.mat'); % get all the .mat files in the current folder
files = files(1:end-2);
load actiCap64; 
protocolType = 'SF_ORI';
projectName = 'ADGammaProject';
outputFolder = 'D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10\interpolatedData\decimatedData\';
%%
for i = 1:length(files)
    load(files(i).name); % load the current file
    
    badElecsTotal = unique([badElecs.badImpedanceElecs;badElecs.noisyElecs;badElecs.flatPSDElecs]);
    goodElecs = setdiff(1:64, badElecsTotal);    
    
    tmp = eeg_emptyset();
    tmp.data = eegData;
    tmp.chanlocs = chanlocs;
    tmp.srate = 250;
    tmp = eeg_checkset(tmp);
    tmp = pop_select(tmp,'nochannel',badElecsTotal);
    % Create the spherical interpolation weights
    if size(eegData,1) == 64
    interpWeights = pop_interp(tmp,chanlocs,'spherical');
    eegData = interpWeights.data;

%     for j = 1:length(badElecs)
%         eegDataInterp(badElecs(j),:,:) = interpWeights(j,:) * eegData(goodElecs,:,:);
%     end
    
    % Save the interpolated data to a new file
%     save(strcat(files(i).name(1:end-4),'_ip.mat'), 'eegDataInterp');
%     save( fullfile(outputFolder, 'eegDataInterp');
    save(fullfile(outputFolder,projectName,protocolType,strcat(files(i).name(1:end-4),'.mat')), 'eegData','badElecs','trialConditionVals','timeVals');
   
    end
end

