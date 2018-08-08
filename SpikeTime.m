function [SpkTime] = SpikeTime(StimOrder,SpikeData,nRep,Time,Win);
%Function to extract spike times separated by frequency and amplitude for
%pure tone stimuli.
%
%Inputs:
%   StimOrder = 3xN array containing tone onset time (in seconds), tone frequency,
%               and tone amplitude
%   SpikeData = matrix containing [Time Stamp; Time from start of experiment;
%                                   Time from start of trial; Trial #]
%   nRep = number of stimulus repeats
%   Win = Size of window from onset of tone for which to look for spikes
%         (in seconds)
%
%Outputs:
%   SpkTime = 3-D cell array (nFreq x nAmp x nRep) containing separated
%             spike times.


%Find number of spikes per trial in time window of tone onset
ResponseT = cell(nRep,length(Time));
for i = 1:nRep
    a = find(SpikeData(4,:) == i); %Find spikes that occur in trial i
    SpikeT = SpikeData(3,a); %Time of spikes that occur within trial j
    for j = 1:length(Time)
        b = find(SpikeT > (Time(j) - Win) & SpikeT <= (Time(j) + Win));
        if ~isempty(b)
            ResponseT{i,j} = SpikeT(b) - Time(j); %Time of spike occurance relative to tone onset.
        end
    end   
end 


%Separate spike times by tone frequency and amplitdue
fList = sort(unique(StimOrder(2,:)));
aList = sort(unique(StimOrder(3,:)));
SpkTime = cell(length(fList),length(aList),nRep);
for i = 1:length(StimOrder)
    a = find(fList == StimOrder(2,i));
    b = find(aList == StimOrder(3,i));
    for j = 1:nRep    
        SpkTime{a,b,j} = [SpkTime{a,b,j} ResponseT{j,i}];
    end
end