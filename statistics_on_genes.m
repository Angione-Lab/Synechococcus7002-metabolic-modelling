load('SynechococcusPCC7002.mat'); %fbamodel
load('Syn7002_IDs.mat'); % list of gene IDs extracted from transcriptomic reads file

genes = fbamodel.genes; 
genes_in_dataset = Syn7002_IDs;


V = numel(genes);
M = 2; %number of objectives


%% COMPUTE PROFILES OF THE NONDOMINATED POINTS
%n_population = non_dominated(:,M+2);
%n_individual = non_dominated(:,M+3);
%
% profiles_nondom = zeros(numel(n_population),V);
% for i = 1:numel(n_population)
%     file_name = ['populations_recon2_merged_bio_IDH1-2/population' num2str(n_population(i)) '.mat'];
%     file_pop = load(file_name);
%     chromosome = file_pop.chromosome;
%     profiles_nondom(i,:) = chromosome(n_individual(i),1:V);
% end
%load('non_dominated_chromosomes.mat'); %if the profiles_nondom have been already computed before
%load('others_chromosomes.mat');
load('transcripts.mat')
all_objpairs = transcripts'; % choose dataset
all_solutions = transcripts'; % same dataset here

%%
all_biomass_values = all_solutions(:,1);

% lim_inf = mean(all_biomass_values) - 2*std(all_biomass_values);
% lim_inf = 1e-10;
% lim_sup = mean(all_biomass_values) + 2*std(all_biomass_values);
% 
% ix_low_biomass = find(all_biomass_values < lim_inf);
% ix_high_biomass = find(all_biomass_values > lim_sup);
% 
% ix_of_interest = ix_high_biomass; disp('I''m doing statistics on the high-biomass points');
% ix_of_interest = ix_low_biomass; disp('I''m doing statistics on the low-biomass points');

profiles = all_objpairs;%(ix_of_interest,:);
genes_vs_profiles = profiles'; %we need to use profiles' and not profiles, because otherwise it would compute the correlation (and all the following measures) between profiles along all the genes, while we want the correlation between genes along the profiles


dist_correlation_vector = pdist(zscore(genes_vs_profiles), 'correlation'); %pdist_corr_1 is my modified version of pdist, in which we use x - 1 instead of x - mean(x) all over the place in the definition of the correlation distance here  http://www.mathworks.co.uk/help/stats/pdist.html?refresh=true
dist_correlation_matrix = squareform(dist_correlation_vector);

clusterTree = linkage(dist_correlation_vector, 'average');







%% *******************************************************************************************************************************************************************************
%try k-means clustering with a few number of clusters: the result was
%that the max(silhouette) was k=10 for the high_biomass_points and k=6 for the low_biomass_points
%Note that we discarded trivial clustering solutions with too few clusters

% prompt = 'k-means: Press ''y'' if the number of cluster is known, or any other key to execute silhouette analysis ';
% answer = input(prompt,'s');
% 
% if strcmp(answer,'y')
%     mean_silhouette = zeros(1,30);
%     for NoClust = 2:30
%         [cidx, ctrs] = kmeans(genes_vs_profiles, NoClust, 'dist','correlation', 'rep', 5, 'disp','final'); %kmeans_corr_1 is my modified version of pdist, in which we use x - 1 instead of x - mean(x) all over the place in the definition of the correlation distance here  http://www.mathworks.co.uk/help/stats/pdist.html?refresh=true
%         %figure;
%         
%         % We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
%         % The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are.
%         figure;
%         [silh5,h] = silhouette(genes_vs_profiles,cidx,'corr');
%         mean_silhouette(NoClust) = mean(silh5);
%         h = gca;
%         h.Children.EdgeColor = [.8 .8 1];
%         close(gcf) %we close the silhouette plot otherwise we have one plot per each index of the loop
%         xlabel 'Silhouette Value';
%         ylabel 'Cluster';
%         
%     end
% end
% 
% figure
% plot(mean_silhouette);
%%
 prompt = 'k-means: what is the number of clusters chosen after inspection of the mean_silhouette plot? (10 for high biomass, 9 for low biomass) ';
 num_clusters = input(prompt)


%% k-means clustering

[cidx, ctrs] = kmeans(genes_vs_profiles, num_clusters, 'dist','cityblock', 'rep',5, 'disp','final'); %kmeans_corr_1 is my modified version of pdist, in which we use x - 1 instead of x - mean(x) all over the place in the definition of the correlation distance here  http://www.mathworks.co.uk/help/stats/pdist.html?refresh=true
figure
for c = 1:num_clusters
    subplot(5,ceil(num_clusters/5),c);
    plot(genes_vs_profiles((cidx == c),:)');

    %axis tight
end
suptitle('K-Means Clustering of Genes');

% We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
% The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are.
figure;
[silh5,h] = silhouette(genes_vs_profiles,cidx,'cityblock'); %silhouette_corr_1 is a special version of silhouette
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';


% PLOT AVERAGE IN EVERY CLUSTER
figure
for c = 1:num_clusters
    subplot(5,ceil(num_clusters/5),c);
    plot(ctrs(c,:)');
    %axis tight
    axis off    % turn off the axis
end
suptitle('K-Means Clustering of Profiles');
%
%clustergram(genes_vs_profiles(:,2:end));


%% *******************************************************************************************************************************************************************************
% try hierarchical clustering with a few number of clusters: the result was
% that the max(silhouette) was with k=12 for the high_biomass_points, while
% the max(silhouette) was with k=6 for the low_biomass_points
% Note that we discarded trivial clustering solutions with too few clusters

%  prompt = 'Hierarchical clustering: Press ''y'' if the number of cluster is known, or any other key to execute silhouette analysis ';
%  answer = input(prompt,'s');
% 
% if strcmp(answer,'y')
%     mean_silhouette = zeros(1,30);
%     for NoClust = 2:30
%         clusters = cluster(clusterTree, 'maxclust', NoClust);
%          
%         %We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
%         %The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are.
%         figure;
%         [silh5,h] = silhouette(genes_vs_profiles,clusters,'corr'); %silhouette_corr_1 is a modified version
%         mean_silhouette(NoClust) = mean(silh5);
%         h = gca;
%         h.Children.EdgeColor = [.8 .8 1];
%         close(gcf) %we close the silhouette plot otherwise we have one plot per each index of the loop
%         xlabel 'Silhouette Value';
%         ylabel 'Cluster';
%         
%     end
%     figure
%     plot(mean_silhouette);
% end
% 
% prompt = 'Hierarchical clustering: what is the number of clusters chosen after inspection of the mean_silhouette plot? (12 for high biomass, 6 for low_biomass)';
% num_clusters = input(prompt)
% 
% %% Hierarchical clustering
% % 
% clusters = cluster(clusterTree, 'maxclust', num_clusters);
% numel_clusters = accumarray(clusters,1);   %cardinality of clusters, i.e. number of genes (elements) in each cluster
% gene_names_in_cluster = cell(num_clusters,1);

% figure
% for c = 1:num_clusters
%     subplot(3,ceil(num_clusters/3),c);
%     plot(genes_vs_profiles((clusters == c),:)');
%     genes_in_cluster = find(clusters == c);
%     gene_names_in_cluster{c} = strjoin(fbamodel.genes(genes_in_cluster)',', ');
%     title({   [num2str(numel_clusters(c)) ':  ' strjoin(fbamodel.genes(genes_in_cluster(1:min(numel_clusters(c),3)))',', ')]     ,     strjoin(fbamodel.genes(genes_in_cluster(4:min(numel_clusters(c),6)))',', ')     }, 'FontSize',7);
%     %axis tight
% end
% suptitle('Hierarchical Clustering of Genes');
% 
% % We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
% % The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are.
% figure;
% [silh5,h] = silhouette(genes_vs_profiles,clusters,'corr');
% h = gca;
% h.Children.EdgeColor = [.8 .8 1];
% close(gcf)
% xlabel 'Silhouette Value';
% ylabel 'Cluster';

% 
% 
% % 
% % %% PLOT DISTANCE BETWEEN GENES THROUGH MULTI-DIMENSIONAL SCALING APPLIED TO THE CORRELATION MATRIX TO MAKE IT BECOME A DISTANCE MATRIX
% % 
%[Y,stress] = mdscale_robust(dist_correlation_vector,2,'sstress','metricsstress');
%[Y,stress] = mdscale_robust(dist_correlation_vector,2,'criterion','metricstress');
options = statset('MaxIter',500);
[Y,stress] = mdscale_robust(dist_correlation_vector,2,'criterion','sstress','start','random','Options',options); %mdscale_robust is a variation which circumvents colocation error by multiplying dissimilarities by a scalar

% figure
% plot(Y(:,1),Y(:,2),'.','LineWidth',2);
% C = clusters; %colour according to hierarchical clustering
% colormap(jet(256))
% scatter(Y(:,1),Y(:,2),200,C,'.');
% title(['Hierarchical clustering (k=' num2str(numel(clusters)) ')']);
% labels = num2str((1:size(Y,1))','%d');    %'
% text(Y(:,1), Y(:,2), labels, 'horizontal','left', 'vertical','bottom')
% colorbar;
% gname(genes);

figure
C = cidx;  %colour according to kmeans clustering
colormap(jet(256))
scatter(Y(:,1),Y(:,2),200,C,'.');
title(['K-Means Clustering (k=' num2str(numel(unique(cidx))) ')']);
labels = num2str((1:size(Y,1))','%d');    %'
text(Y(:,1), Y(:,2), labels, 'horizontal','left', 'vertical','bottom')

%colorbar;
