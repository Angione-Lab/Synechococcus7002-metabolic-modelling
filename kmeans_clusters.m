%% k-means clustering of the influenza flux data for each virus strain
%load('D:\E Documents\Influenza\Clusters\kmeans\H3N2.mat');
%load('D:\E Documents\Influenza\Clusters\kmeans\H5N1.mat');
%load('D:\E Documents\Influenza\Clusters\kmeans\H7N7.mat');
%load('D:\E Documents\Influenza\Clusters\kmeans\H7N9.mat');
%load('D:\E Documents\Influenza\Clusters\kmeans_all\VirusFlux.mat');
load 'all_data_kmeans_14.06.mat'
% A struct containing the flux values for each virus, rows are fluxes
% columns are time samples (7785 x 4) observations x variables


%% To help process large dataset, invoke a parallel pool of workers. Specify options for paralell
% computing. 
% pool = parpool;                      % Invokes workers
% stream = RandStream('mlfg6331_64');  % Random number stream
% options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);

%% initialise variables
% set number of clusters and observations 
numClust = 5;
numObs = 24;

% set fluxes to use
%transcripts % use transcripts values
%fluxes = VirusFlux; % get flux values

% create an array to store the total sum of squared errors
%VirusFlux_SSE = zeros(numClust, 1);
transcripts_SSE = zeros(numClust, 1);
% create an array to store cluster allocations
%VirusFlux_kmeans_clusters = zeros(numObs,numClust);
transcripts_kmeans_clusters = zeros(numObs,numClust);
% Create a new array to hold the silhouette values
silh_values = zeros(numObs,numClust); %rows = fluxes

%% Run kmeans using defualt distance (sqeuclidean)
% fluxes/transcripts is the data matrix 
% i is the number of clusters to try
% options is set above for parallel computing
% MaxIter - maximum number of iterations
% Display - displays the results of the final iteration
% Replicates - number of times to repeat clustering using the new initial
% cluster centroid positions
% idx is a column vector showing the cluster number that each object has
% been assigned to
% C contains the cluster centroid locations
% sumd contains within-cluster sums of point-to-centroid distances
% D contains distances from each point to every centroid

for i=2:numClust
    [idx,C,sumd,D] = kmeans(transcripts,i,'MaxIter',100,'Display','off ','Replicates',10);
    s = silhouette(transcripts,idx);        
    silh_values(:,i) = s;

    % store cluster allocations
    %VirusFlux_kmeans_clusters(:,i) = idx;
    transcripts_kmeans_clusters(:,i) = idx;
    %store sum of squared error 
    %VirusFlux_SSE(i) = sum(sumd);
    transcripts_SSE(i) = sum(sumd);
end
            

% Store the mean Silhouette value. Values close to 1 are desirable as it indicates that points are very distant from neighbouring clusters.
%VirusFlux_avgSilh_kmeans = nanmean(silh_values); % average sihouette value for each cluster      
transcripts_avgSilh_kmeans = nanmean(silh_values); % average sihouette value for each cluster      

% Plot the mean silhouette
figure
plot(transcripts_avgSilh_kmeans);
      
% Plot clusters?



     
 



