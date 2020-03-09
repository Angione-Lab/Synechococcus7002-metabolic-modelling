% LASSO regularization with alpha = 1 which returns fitted least-squares regression coefficients for linear models of
% transcript/flux/transcript+flux data (x) and the growth rates (y) in 12 growth conditions

load('transcripts_subset.mat') % load transcript data corresponding to 12 growth conditions
load('all_atp_flux_subset.mat') % load flux data corresponding to 12 growth conditions
load('all_p1_flux_subset.mat')
load('all_p2_flux_subset.mat')
load('ATPTF_subset.mat') % load transcript + flux data corresponding to 12 growth conditions
load('P1TF_subset.mat')
load('P2TF_subset.mat')
load('Y2.mat') % load growth rates corresponding to 12 growth conditions

[B_transcripts,fitInfo_transcripts]=lasso(transcripts_subset,Y2);

[B_ATP,fitInfo_ATP]=lasso(all_atp_flux_subset,Y2);
[B_P1,fitInfo_P1]=lasso(all_p1_flux_subset,Y2);
[B_P2,fitInfo_P2]=lasso(all_p2_flux_subset,Y2);

[B_ATPTF,fitInfo_ATPTF]=lasso(ATPTF_subset,Y2);
[B_P1TF,fitInfo_P1TF]=lasso(P1TF_subset,Y2);
[B_P2TF,fitInfo_P2TF]=lasso(P2TF_subset,Y2);
