function [STRFclusterdata] = calcSTRFcluster(STRF,randSTRF,sigma_b,p,tC)
%Identify lobes in STRF using Stat4ci toolbox
%
%INPUT:
%   STRF = unsmoothed STA array
%   randSTRF = random unsmoothed STA to serve as noise
%   sigma_b = smoothing window
%   p = significance value
%   tC = cluster threshold
%
%OUTPUT:
%   STRFclusterdata = significant clusters found

%Calculate z-scored STRF
smoothSTRF = imgaussfilt(STRF,sigma_b);
smoothrandSTRF = imgaussfilt(randSTRF,sigma_b);
Zsta = ZScoreSCi(smoothSTRF,[mean(smoothrandSTRF(:)), std(smoothrandSTRF(:))]);

wholePic = ones(size(Zsta,1),size(Zsta,2));
clusters = StatThresh(Zsta, p, sigma_b, tC, wholePic);

%Identify significant clusters
[STRFclusterdata.tCi, STRFclusterdata.info] = DiplayRes_STRF(clusters,smoothSTRF);
