function [mask, clustmean] = matchSTRFclust(STRFmaskREF,STRFmaskCOMP, strfCOMP)
%Generates new mask for second input to match the cluster labels to the
%first input.
%
%Inputs:n (all from calcSTRFcluster.m)
%   STRFmaskREF = cluster mask with unique integer labels for each cluster
%   STRFmaskCOMP = cluster mask with unique integer labels for each cluster
%   strfCOMP = COMP STRF (optional input - calc mean COMP STRF value using REF mask)
%
%Outputs:
%   mask = new cluster mask for STRFmaskCOMP with unique integer labels for
%          each cluster to match labels in STRFmaskREF if clusters overlap
%   clustmean = [cluster ID, mean cluster using REF mask]

mask = zeros(size(STRFmaskCOMP));
nClustREF = max(unique(STRFmaskREF(:)));
nClustCOMP = max(unique(STRFmaskCOMP(:)));
clustmean = NaN;

clustcnt = 0;
counter = 0;
for j = 1:nClustCOMP  
    tempmaskCOMP = STRFmaskCOMP;
    tempmaskCOMP(tempmaskCOMP ~= j) = 0;        
    tempmaskCOMP(tempmaskCOMP > 0) = 1;
    clustsizeCOMP = sum(tempmaskCOMP(:));
    match = 0;
    for i = 1:nClustREF
        tempmaskREF = STRFmaskREF;
        tempmaskREF(tempmaskREF ~= i) = 0;            
        tempmaskREF(tempmaskREF > 0) = 1;
        clustsizeREF = sum(tempmaskREF(:));
        overlap = sum(sum(tempmaskREF.*tempmaskCOMP));
        if overlap > 0.5*min(clustsizeREF,clustsizeCOMP)
            mask = mask + i.*tempmaskCOMP;
            clustcnt = clustcnt + 1;
            match = 1;
            if nargin > 2 %Mean COMP STRF value using REF mask
                clustmean(clustcnt,1) = i; %cluster ID
                [~,~,val] = find(tempmaskREF.*strfCOMP);
                clustmean(clustcnt,2) = mean(val);
            end
        end
    end

    if match == 0 %Re-label clusters that do not match
        counter = counter+1;
        mask = mask + (counter+nClustREF).*tempmaskCOMP;
    end
end    
    


