%This script calculates Pearson correlation coefficients 
%between a subset of the transcript/flux FC dataset and
%vector of growth rates across 12 conditions

%Create output vectors to store correlation coefficients, p-values, and lower and upper bounds of confidence
%intervals, changing no. of rows for transcripts(3187), fluxes (742), or both (3929)
corr=zeros(3187,1); % PCC
pval=zeros(3187,1); % p-value
lb95=zeros(3187,1); % lower bound for 95% confidence
ub95=zeros(3187,1); % lower bound for 95% confidence

 N = size(transcripts,2); % Specify size and type of data matrix: transcripts/all_atp_flux/all_p1_flux,all_p2_flux
 for i=1:N
     [R,P,RL,RU] = corrcoef(transcripts(:,i),Y2); %Y2 contains growth rates
     corr(i) = R(1,2);
     pval(i) = P(1,2);
     lb95(i) = RL(1,2);
     ub95(i) = RU(1,2);
 end
  
%Calculate coefficient of determination for transcripts/fluxes and growth rates
 
%  for i=1:N
%  rsquared=rsquare(transcripts(:,i),Y2);
%  rsq(i)= rsquared;
%  end
