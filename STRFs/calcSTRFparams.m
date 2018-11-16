function STRFclusterdata = calcSTRFparams(STRFclusterdata,STRF,times, freqs)
%Calculate size of negative and positive lobes (if they exist), peak, etc.
%using results of CLUSTER test
%
%INPUTS:
%   STRFclusterdata = data calculated from calcSTRFcluster.m
%   params = contains frequency and time bins 
%
%OUTPUTS:
%   STRFparams = cluster parameters and parameter labels

if STRFclusterdata.info.Clusterdata(end) > 0 
    
    clustIDs = unique(STRFclusterdata.tCi.ClusMask(:))';
    clustIDs = clustIDs(clustIDs > 0);
    counter = 0;
    for i = clustIDs
        counter = counter + 1;
        %Make binary mask for each cluster
        tempMask = STRFclusterdata.tCi.ClusMask;
        tempMask(tempMask ~= i) = 0;
        ClustMask = tempMask;
        ClustMask(ClustMask > 0) = 1;

        %Calculate parameters for each cluster
        [f, t, x] = find(ClustMask.*STRF);

        deltaf = freqs(max(f))-freqs(min(f));
        deltat = times(max(t))-times(min(t));

        [~, j] = max(x);

        fCenter = freqs(f(j));
        tDelay = times(t(j));
        STRFclusterdata.params.data(counter,:) = [i, deltat, deltaf/1000, tDelay, fCenter/1000, sum(ClustMask(:)), mean(x)];
    end
     STRFclusterdata.params.labels = {'clust ID','time width (s)','freq width (kHz)','peak time (s)','peak freq (kHz)','size (pixels)', 'mean value'};
else 
     STRFclusterdata.params.data = [0 NaN NaN NaN NaN 0 NaN];
     STRFclusterdata.params.labels = {'clust ID','time width (s)','freq width (kHz)','peak time (s)','peak freq (kHz)', 'size (pixels)', 'mean value'};
end