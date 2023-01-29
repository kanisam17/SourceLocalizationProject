% folderLORETA contains all text files in sLORtoText folder.
% 
% 
folderStr = 'sLORETA_Thres10';
folderLORETA = fullfile('N:\Projects\Kanishka_SourceLocalizationProject\data\decimatedData\',folderStr); %change folders
folderOutput = fullfile('D:\OneDrive - Indian Institute of Science\Supratim\Projects\Kanishka_SourceLocalizationProject\data\',folderStr);
clc;clear
% % load('caseListAgeMatched.mat') % for caseControl
folderStr = 'newData';
folderLORETA = fullfile('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\caseControl\',folderStr); %change folders
folderOutput = fullfile('D:\Kanishq\NewProject\TLSAEEGProjectPrograms\decimatedData\LORETA\caseControl\newData',folderStr);
strList = {'control','case'};
saveRawLORETAData(caseList,strList,folderLORETA,folderOutput);



% %% for control averaged against each case
% for i = 1:17
%     filename = sprintf('avgControl%d.mat', i);
%     mDataBL = caseList.subjectNameListMean{1, i}.mDataBL;
%     mDataST = caseList.subjectNameListMean{1, i}.mDataST;
%     pVals = caseList.subjectNameListMean{1, i}.pVals;
%     tStats = caseList.subjectNameListMean{1, i}.tStats;
%     save(filename, 'mDataBL', 'mDataST', 'pVals', 'tStats');
% end