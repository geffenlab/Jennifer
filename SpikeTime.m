function [SpkTime_Laser, SpkTime_NoLaser, fList, aList] = SpikeTime(afo,SpikeData,nRep,dur,Win);
%afo = Information about the order of stimulus presentation
%SpikeData = matrix containing [Time Stamp; Time from start of experiment;
%Time from start of trial; Trial #]
%nRep = number of trials
%Win = Size of window from onset of tone for which to look for spikes 
%dur = stimulus duration (of a single repetition)
Time = 0:0.5:dur; %Each repetition is 400 seconds long, each trial is 500ms long
Time = Time(1:end-1);

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

%Assume the laser is on for every other tone, starting with the second tone
StimOrder_Laser = [Time(2:2:end); afo.freqOrder(2:2:end); afo.ampOrder(2:2:end)]; 
StimOrder_NoLaser = [Time(1:2:end); afo.freqOrder(1:2:end); afo.ampOrder(1:2:end)];
ResponseT_Laser = ResponseT(:,2:2:end);
ResponseT_NoLaser = ResponseT(:,1:2:end);

fList = sort(unique(afo.freqOrder));
aList = sort(unique(afo.ampOrder));

SpkTime_Laser = cell(length(fList),length(aList),nRep);
for i = 1:length(StimOrder_Laser)
    a = find(fList == StimOrder_Laser(2,i));
    b = find(aList == StimOrder_Laser(3,i));
    for j = 1:nRep    
        SpkTime_Laser{a,b,j} = ResponseT_Laser{j,i};
    end
end

SpkTime_NoLaser = cell(length(fList),length(aList),nRep);
for i = 1:length(StimOrder_Laser)
    a = find(fList == StimOrder_NoLaser(2,i));
    b = find(aList == StimOrder_NoLaser(3,i));
    for j = 1:nRep
       SpkTime_NoLaser{a,b,j} = ResponseT_NoLaser{j,i};
    end 
end