 function [SignificantVariables, clusters,clustersize,StatMapPermPV] = permutation_cluster_1sample_alld(data, nperm, cluster_th, significance_th,tail,StatMapPermPV)

% function [SignificantVariables, clusters,clustersize,StatMapPermPV] = permutation_cluster_1sample_alld(data, nperm, cluster_th, significance_th,tail,StatMapPermPV)

%

% Performs one-sided (>0) or two-sided cluster-size test on the 'data'.

% The data array can have any number of dimensions (1D, 2D, 3D, 4D etc) in

% the order [observations x variable1 x variable2 x ...]

% Data from each observation are randomly multiplied by +-1 to create

% permutation samples. These samples are converted to pvalues and then the

% cluster_th threshold (in pvalue units) is applied to identify suprathreshold clusters. The

% distribution of the size of suprathreshold clusters is used to assign statistical

% significance to the clusters (also in pvalue units) of the original data.

%

% INPUT:

%   data: observations x variable1 x variable2 x variable3 x ... (supports all dimemsions)

%   nperm: number of permutations

%   cluster_th: cluster defining threshold (in pvalue units)

%   significance_th: significance threshold (alpha value)

%   tail: string 'right' or 'both'. Default is 'both'.

%   StatMapPermPV: (optional) permutation pvalue map (see output)

%

% % for example:

%   load data;

%   nperm = 1000;

%   cluster_th = 0.05;

%   significance_th = 0.05;

%   tail = 'both';

%   StatMapPermPV: permutation pvalue map (nperm x variable1 x variable2 x...)

%

% OUTPUT:

%   SignificantVariables:  times with significant clusters

%   clusters:           clusters above cluster defining threshold

%   clustersize:        size of identified culsters

%

% Author: Dimitrios Pantazis, December 2015

%



%decide one-sided (right) or two-sided (both) test

if ~exist('tail') || strcmp(tail,'right') %if two tail t-test

    func = '';

else %if two-sided t-test

    func = 'abs';

end



%initialize

N = ndims(data); %data is observations x variable1 x variable2 x ... (supports all dimensions)

nobservations = size(data,1);

for n = 2:N

    nvariable(n-1) = size(data,n);

end

cln = repmat({':'},1,N-1); %to select the N-1 dimensions



%create permutation samples and convert them to pvalues (perms x variable1 x variable2 x ...)

if ~exist('StatMapPermPV') %if pvalues have not been precomputed

    

    

    StatMapPerm = single(zeros([nperm nvariable]));

    %first permutation sample is original data

    StatMapPerm(1,cln{:}) = mean(data,1) ./ std(data);

    

    %perform permutations

    for i = 2:nperm %par

        if ~rem(i,100)

            disp(['Create permutation samples: ' num2str(i) ' out of ' num2str(nperm)]);

        end

        perm = single(sign(rand(nobservations,1)-0.5));

        %data_perm = gmultiply(perm,data); %multiply by filling extra dimesions (data_perm = repmat(perm,1,ntimes1,ntimes2) .* data)

        data_perm = repmat(perm,1,nvariable(1)) .* data; %
        
        %permsample = mean(data_perm,1) ./ std(data_perm);

        %StatMapPerm(i,:) = permsample(:); %does this for multiple dimensions: 
        StatMapPerm(i,cln{:}) = mean(data_perm,1) ./ std(data_perm); 

        

    end

    

    %convert to pvalues

    eval([ 'StatMapPermPV = (nperm+1 - tiedrank(' func '(StatMapPerm)))/nperm;' ]);

    

end

clear StatMapPerm;



%find maximum cluster for all permutation samples

[clustermax(1),nclusters,clusters,clustersize] = find_clusters_alld(squeeze(StatMapPermPV(1,cln{:})<=cluster_th) , cluster_th);

for i = 2:nperm

    if ~rem(i,100)

        disp(['Compute max cluster: ' num2str(i) ' out of ' num2str(nperm)]);

    end

    [clustermax(i)] = find_clusters_alld(StatMapPermPV(i,cln{:})<cluster_th,cluster_th);

end



%find cluster threshold

clustermax_sorted = sort(clustermax);

th = clustermax_sorted(nperm*(1-significance_th));



%find significant variables

% SignificantTimes = [clusters{clustersize>th}];

if length(nvariable) == 1

    SignificantVariables = zeros(nvariable,1);

else

    SignificantVariables = zeros(nvariable);

end

SignificantVariables([clusters{clustersize>th}]) = 1;

%figure;imagesc(SignificantVariables);set(gca,'YDir','normal');



%clear pool



