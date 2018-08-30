function [mask, something] = matchSTRFclust(STRFmaskREF,STRFmaskCOMP, ZSCiCOMP)
%Generates new mask for second input to match the cluster labels to the
%first input.
%
%Inputs:
%   STRFmaskREF = cluster mask with unique integer labels for each cluster
%   STRFmaskCOMP = cluster mask with unique integer labels for each cluster
%
%Outputs:
%   mask = new cluster mask for STRFmaskCOMP with unique integer labels for
%          each cluster to match labels in STRFmaskREF if clusters overlap

mask = zeros(size(STRFmaskCOMP));
nClustREF = max(unique(STRFmaskREF(:)));
nClustCOMP = max(unique(STRFmaskCOMP(:)));

if nClustREF >= nClustCOMP
    counter = 0;
    for i = 1:nClustREF
        STRFmaskREF(STRFmaskREF ~= i) = 0;
        tempmaskREF = STRFmaskREF;
        tempmaskREF(tempmaskREF > 0) = 1;
        clustsizeREF = sum(tempmaskREF(:));
        for j = 1:nClustCOMP
            STRFmaskCOMP(STRFmaskCOMP ~= j) = 0;
            tempmaskCOMP = STRFmaskCOMP;
            tempmaskCOMP(tempmaskCOMP > 0) = 1;
            clustsizeCOMP = sum(tempmaskCOMP(:));
            overlap = sum(sum(tempmaskREF.*tempmaskCOMP));
            if overlap > 0.5*min(clustsizeREF,clustsizeCOMP)
               mask = mask + i.*tempmaskCOMP;
            else %If none of the COMP clusters match to a REF cluster
                counter = counter + 1;
                mask = mask + (counter*nClustREF).*tempmaskCOMP;
            end

        end

    end
else
    counter = 0;
    for j = 1:nClustCOMP       
        STRFmaskCOMP(STRFmaskCOMP ~= j) = 0;
        tempmaskCOMP = STRFmaskCOMP;
        tempmaskCOMP(tempmaskCOMP > 0) = 1;
        clustsizeCOMP = sum(tempmaskCOMP(:));
        
        for i = 1:nClustREF
            STRFmaskREF(STRFmaskREF ~= i) = 0;
            tempmaskREF = STRFmaskREF;
            tempmaskREF(tempmaskREF > 0) = 1;
            clustsizeREF = sum(tempmaskREF(:));
            overlap = sum(sum(tempmaskREF.*tempmaskCOMP));
            if overlap > 0.5*min(clustsizeREF,clustsizeCOMP)
               mask = mask + i.*tempmaskCOMP;
            else %Re-label clusters that do not match
                counter = counter+1;
                mask = mask + (counter+nClustREF).*tempmaskCOMP;
            end
        end
    end    
    
end

%Calculate mean activity for ref clusters in comp
if nargin > 2
    
end

