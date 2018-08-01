function STA = genSTRF(spikes,StimParams,win,fig)

if StimParams.chordDuration > 0.005
    rin = StimParams.chordDuration;
    rout = 0.005;
    STIM = cleanResample(StimParams.dbs,rin,rout);
    fps = 1/rout;
else STIM = StimParams.dbs;
    fps = 1/StimParams.chordDuration;
end

%Input cannot include first win ms. 
%ind = find(spikes < win);
%spikesSTA = spikes;
%spikesSTA(ind) = [];

[STA, ~] = genSTA(spikes,STIM,win,fps);

if nargin > 3
    t = [0:20]; %Not sure what these time bins should be...
    plotSTA(StimParams.freqs,t,STA);
end
