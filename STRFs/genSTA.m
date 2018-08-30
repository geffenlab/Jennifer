function [sta, nSpikes] = genSTA(spikes,S,w,fps);
% This function generates a spike-triggered average
%
%Inputs:
%   spikes = spike times (needs to be in seconds from the start of
%            the stimulus
%   S = stimulus spectrogram of m frequencies by n time bins (dbs matrix)
%   w = length of the stimulus window to consider for the STA in seconds (how
%       far back in time do you want your STA to go) (100 ms)
%   fps = stimulus frame rate (1 / chord length)

% convert spike times to stim bins
spikeT = ceil(spikes*fps);

sta = zeros(size(S,1),w/(1/fps)+1);
for i = 1:length(spikeT)
    spikeStim = S(:,spikeT(i) - ((w*fps)):spikeT(i));
    sta = sta + spikeStim;
end

nSpikes = length(spikeT);
sta = sta / nSpikes;