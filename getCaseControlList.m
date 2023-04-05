function subjectNameListFinal = getCaseControlList (subjectList64Chan,cdrList,ageList,genderList)
% function made to create case control list subjects with 64 channels only
% exclusive for source localization project.
%% Add line below in runDisplayData code
%%% subjectList64Chan = uniqueSubjectNames0(goodIndices);
%%
ageLim = 1;
healthyPos = strcmp(cdrList,'HV');
casePos = ~healthyPos;
caseList = setdiff(subjectList64Chan(casePos),[{'217SK'} {'225SK'}]);
controlList = []; controlListCaseNumber = [];
for i=1:length(caseList)
    subjectName = caseList{i};
    pos = find(strcmp(subjectName,subjectList64Chan));
    age = ageList(pos); gender = genderList(pos);
    ageMatchPos = (ageList<= age+ageLim) & (ageList>= age-ageLim);
    genderMatchPos = strcmp(gender,genderList);
    controlPos = healthyPos & ageMatchPos & genderMatchPos;
    controls = subjectList64Chan(controlPos);
    controlList = cat(2,controlList,controls);
    controlListCaseNumber = cat(2,controlListCaseNumber,i+zeros(1,length(controls)));
    disp([num2str(i) '. ' subjectName ' (' num2str(age) ',' gender{1} '): ' num2str(length(controls)) ' controls.']);
end
% % Method 1 for unmatched age/gender
% subjectNameList{1} = unique(controlList); strList{1} = 'Controls';
% subjectNameList{2} = caseList; strList{2} = 'Cases';

% Method 2 for matched age/gender
casesWithControls = unique(controlListCaseNumber);
numValidCases = length(casesWithControls);
subjectNameListFinal = cell(1,2);
for i=1:numValidCases
    subjectNameListFinal{1}{i} = controlList(casesWithControls(i)==controlListCaseNumber);
    subjectNameListFinal{2}{i} = caseList(casesWithControls(i));
end
strList{1} = 'control';
strList{2} = 'case';

end
