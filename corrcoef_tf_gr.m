% This script calculates Pearson correlation coefficients 
% between a subset of the transcript/flux FC dataset and
% vector of growth rates across 12 conditions
% N = size(minNormATPTF,2); % Specify size of data matrix
% for i=1:N
%     tmp = corr(minNormATPTF(:,i),Y1,'Type','Pearson'); %Y1 contains doubling %, Y2 contains growth rates
%     Ct(i,:)=tmp;
% end
    
N = size(absolute_P2TF_FC,2); % Specify size of data matrix
for i=1:N
    tmp = corr(absolute_P2TF_FC(:,i),Y2,'Type','Pearson'); %Y1 contains doubling %, Y2 contains growth rates
    D(i,:)=tmp;
end