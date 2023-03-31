% folderLORETA contains all text files in sLORtoText folder.
% 
% load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\subjectNameListFinal.mat', 'subjectNameListFinal') % load subject list for project 
% load('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\subjectNamesWith64elecs.mat') % load uniqueSubjectList with 64 electrodes only.
folderStr = 'text';
folderLORETA = fullfile('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10\interpolatedData\decimatedData\AdGammaProject\',folderStr); %change folders
folderOutput = fullfile('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\sLORETA_Thres10\interpolatedData\decimatedData\sourceData\LORETA\data\',folderStr);
strList = {'mid','old'};

% % % load('caseListAgeMatched.mat') % for caseControl
% folderStr = 'newData';
% folderLORETA = fullfile('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\caseControl\',folderStr); %change folders
% folderOutput = fullfile('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\caseControl\newData',folderStr);
% strList = {'control','case'};

% %% Check positions and trial index of each subject 
% pos = cell(1, 2); trialInfoPos = cell(1, 2);
% for i = 1:2
%     pos{i} = [];
%     for j = 1:length(subjectNameListFinal{i})
%         filename = subjectNameListFinal{i}{j};
%         idx = find(ismember(subjectNames, filename));
%         if ~isempty(idx)
%             pos{i} = [pos{i} idx];
%         end
%     end
%     trialInfoPos{i} = trialIdxAllSub(pos{i});
% end


%% save files
% saveRawLORETAData(subjectNameListFinal,strList,folderLORETA,folderOutput,trialIdxListFinal);
saveRawLORETAData(subjectNameListFinal,strList,folderLORETA,folderOutput);
clc;clear

% %% for control averaged against each case
% for i = 1:17
%     filename = sprintf('avgControl%d.mat', i);
%     mDataBL = caseList.subjectNameListMean{1, i}.mDataBL;
%     mDataST = caseList.subjectNameListMean{1, i}.mDataST;
%     pVals = caseList.subjectNameListMean{1, i}.pVals;
%     tStats = caseList.subjectNameListMean{1, i}.tStats;
%     save(filename, 'mDataBL', 'mDataST', 'pVals', 'tStats');
% end