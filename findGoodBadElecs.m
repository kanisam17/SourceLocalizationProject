function badElecsPerSubject = findGoodBadElecs(subjectNameListFinal,dataStr,folderSourceString)

x = load('goodProtFlag');

for i=1:2
    subjectNames = subjectNameListFinal{i};
    
    for j=1:length(subjectNames)
        subjectName = subjectNames{j};
        
        disp([ dataStr{i} ', ' subjectName]);
        
        [expDates,protocolNames] = getProtocolDetailsForAnalysis('ADGammaProject',subjectName,'SF_ORI');
        
        pos = strcmp(x.uniqueSubjectNames,subjectName);
        goodProts = x.goodProtFlagList{pos};
        
        for k=1:length(protocolNames)
            if goodProts(k) % If this is a good protocol
                
                
                
            end
        end
    end
end