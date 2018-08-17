function STRFclusterdata = calcSTRFparams(STRFclusterdata,times, freqs)
%Calculate size of negative and positive lobes (if they exist), peak, etc.
%using results of CLUSTER test
%
%INPUTS:
%   STRFclusterdata = data calculated from calcSTRFcluster.m 
%   params = contains frequency and time bins 
%
%OUTPUTS:
%   STRFparams = cluster parameters and parameter labels

[f, t, x] = find(STRFclusterdata.tCi.ClusMask.*STRFclusterdata.tCi.ZSCi);

deltaf = freqs(max(f))-freqs(min(f));
deltat = times(max(t))-times(min(t));

[~, j] = max(x);

fCenter = freqs(f(j));
tDelay = times(t(j));

STRFclusterdata.params.data = [deltat, deltaf/1000, tDelay, fCenter/1000];
STRFclusterdata.params.labels = {'time width (s)','freq width (kHz)','peak time (s)','peak freq (kHz)'};