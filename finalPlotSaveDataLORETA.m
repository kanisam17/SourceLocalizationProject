%% Convert final plotting data to text file: a precursor for LORETA gui to
% make figures

grp = {'mid','old'}; freq = {'A','SG','FG'};

% For each Subjects 
for j = 1:size(grp,2)
    for i = 1:size(freq,2)
        dataLORDeltaP{j}{i} = squeeze(dataDeltaP{j}(i,:,:));
        writematrix(dataLORDeltaP{j}{i}',fullfile(pwd,['dataLORDeltaP' '-' grp{j} '-' freq{i} '.txt']), 'Delimiter' , 'space');
    end
end

%Averaged across subjects
for j = 1:size(grp,2)
    for i = 1:size(freq,2)
        dataLORDeltaP{j}{i} = squeeze(dataDeltaP{j}(i,:,:));
        avgDataLORDeltaP{j}{i} = mean(dataLORDeltaP{j}{i},2);
%         writematrix(dataLORDeltaP{j}{i}',fullfile(pwd,['dataLORDeltaP' '-' grp{j} '-' freq{i} '.txt']), 'Delimiter' , 'space');
        writematrix(avgDataLORDeltaP{j}{i}',fullfile(pwd,['avgDataLORDeltaP' '-' grp{j} '-' freq{i} '.txt']), 'Delimiter' , 'space');
    end
end
