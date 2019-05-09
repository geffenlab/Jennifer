% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear all
seed = round(rand(1)*1000);
rng(seed);
disp(seed)

filename = 'TC004_180816_blueAcuteSpeaker';
fs = 200000;
f = [3000 70000];
stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
minF = f(1); % minimum frequency to test
maxF = f(2); % max frequency to test
octaveSteps = 1/6; % distance between the tones in octaves
nLogSteps = 50; % number of log steps
tDur = 50; % tone duration in ms
ITI = 450; % inter tone interval duration in ms
attenuations = [0];
repeats = 5; % number of repeats of each tone


% Load the calibration filter
filtName = '180816_blueAcuteSpeaker_NIDAQ_3k-80k_fs200k.mat';
load(filtName);

% Create stimuli
totalDur = (ITI+tDur)/1000; % total duration in s
toneDur = tDur/1000;

if length(f)==2
    if strcmp(stepType,'octaves')
        freqs = minF;
        while freqs(end)<maxF
            freqs(length(freqs)+1) = freqs(length(freqs))+(freqs(length(freqs))*octaveSteps);
        end
        freqs = round(freqs(1:end-1));
    elseif strcmp(stepType,'log')
        freqs=exp(linspace(log(minF),log(maxF),nLogSteps));
        freqs = round(freqs);
    end
else
    freqs = f;
end

% preallocate variable
stimArray = zeros(length(freqs)*length(attenuations),round(totalDur*fs));
%events = stimArray;
loc = 1; % placement of the tone in the zeros
ind = 1; % initiate index
for ii = 1:length(freqs)
    for jj = 1:length(attenuations)
        t = generateTone(freqs(ii),(3*pi)/2,toneDur,fs); % Make tone
        t = envelopeKCW(t,5,fs); % envelope
        t = t.*10^(-attenuations(jj)/20); % attenuate
        stimArray(ind,loc:loc+length(t)-1) = t;
%         stimArray(ind,:) = conv(stimArray(ind,:),FILT,'same');
        %events(ind,loc:loc+length(t)-1) = ones(1,length(t))*5;
        index(ind,1) = freqs(ii);
        index(ind,2) = attenuations(jj);
        ind = ind+1;
    end
end

stim = [];
ind = 1; % initiate index
toneOrder=[];

for ii = 1:repeats
    disp(ii)
    reg_ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
    flp_ro = fliplr(reg_ro);
    tmp_allorder_ro = [reg_ro; flp_ro];
    ro = reshape(tmp_allorder_ro,1,size(tmp_allorder_ro,1)*size(tmp_allorder_ro,2));
    
    for idx = 1:length(ro)
    
    %stimT_reg = reshape(stimArray(reg_ro,:)',1,length(stimArray(:)));
    %stimT_flp = reshape(stimArray(flp_ro,:)',1,length(stimArray(:)));
    %stimT_all = [stimT_reg; stimT_flp];
    %stimT = reshape(stimT_all,1,size(stimT_all,1)*size(stimT_all,2));

    %stim((ii-1)*length(stimArray(:))+1:ii*length(stimArray(:))*2,1)=stimT;
    stim = [stim stimArray(ro(:,idx),:)];
    end
    toneOrder = [toneOrder,ro];
end

stim = conv(stim/10,FILT,'same');

% Make stim info
stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.frequencies = f;
stimInfo.stepType = stepType; % log or octaves, if length(f)>3 then those frequencies will be used
stimInfo.octaveSteps = octaveSteps; % distance between the tones in octaves
stimInfo.nLogSteps = nLogSteps; % number of log steps
stimInfo.tDur = tDur; % tone duration in ms
stimInfo.ITI = ITI; % inter tone interval duration in ms
stimInfo.attenuations = attenuations;
stimInfo.repeats = repeats; % number of repeats of each tone
stimInfo.filterName = filtName;
stimInfo.index = index;
stimInfo.order = toneOrder;
stimInfo.stimGenFunc = 'toneSequenceGen.m';
% order = toneOrder;
disp(length(stim)/fs)

%% Write file and save WITHOUT laser and event triggers


cd('C:\Users\Jennifer\Documents\MATLAB\Stimuli\Tuning Curve')
fn = [filename '.wav'];
stim = (stim);
audiowrite(fn, stim, fs)
save([filename '_stimInfo.mat'],'stimInfo')

%% Write file and save WITH laser and event triggers - for blue booth

laserdur = 250; %laser duration
laserON = 508; %400 480 508 %laser onset within 2 clicks 
evpulse = 5;

laser1 = zeros(2*(tDur+ITI)*fs/1000,1);
ev1 = zeros(2*(tDur+ITI)*fs/1000,1);

ev1(1:evpulse*fs/1000) = 0.5;
ev1((tDur+ITI)*fs/1000:(tDur+ITI+evpulse)*fs/1000) = 0.5;
laser1(laserON*fs/1000:(laserON+laserdur)*fs/1000) = 0.5;

nTones = length(stim)/((tDur+ITI)*fs/1000);
ev = repmat(ev1,nTones/2,1);
laser = repmat(laser1,nTones/2,1);

stimInfo.laserdur = laserdur;
stimInfo.laserON = laserON;

cd('C:\Users\Jennifer\Documents\MATLAB\Stimuli\Tuning Curve')
fn = [filename 'C'];
STIM = [stim' ev laser];
audiowrite([fn '.wav'], STIM, fs)
save([fn '_stimInfo.mat'],'stimInfo')





