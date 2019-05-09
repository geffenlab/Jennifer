function SpikeCount = getSpikeCount_Trial(StimOrder, SpikeData, Win)
%Win = vector containing (1) how far back (non-positive value) and (2) how far
%      forward (non-negative value)

Time = StimOrder(1,:);
nTrial = length(Time);
SpikeCount = NaN(nTrial,1);
SpikeT = SpikeData(2,:);
for j = 1:nTrial
    b = find(SpikeT > (Time(j) + Win(1)) & SpikeT <= (Time(j) + Win(2)));
    SpikeCount(j) = length(b);
end
       
        
        