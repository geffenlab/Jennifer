function smoothedFR = smoothFRx3(SpikeData,nRep,binsize,Win,gwin)
%Smooths spike data by convolving firing rates with a gaussian filter
%Input:
    %SpikeData = matrix [time(relative to trial onset); trial #]
    %binsize = size of bin (in seconds) for histogram of spike data
    %gwin = [N,ALPHA] creates a gaussian smoothing window of N points where
    %       ALPHA is the reciprocal of the standard deviation.


t = [Win(1):binsize:Win(2)];

f = histc(SpikeData(1, :), t)/(nRep*binsize); %dividing by nRep averages the number of spikes over the trials instead of just summing them and dividing by binsize makes it spikes per second instead of spikes per 5ms.

t1 = t(1:end-1)+(t(2)-t(1))/2;
f1 = f(1:end-1);

GaussFilter = gausswin(gwin(1), gwin(2)); % 4.8 gives a 5ms half-width.
GaussFilter = GaussFilter/sum(GaussFilter); %Definitely need this!
smoothedFR = conv(f1, GaussFilter, 'same');

end