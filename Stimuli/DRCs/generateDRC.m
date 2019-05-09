%% SETUP
% output directory
cd 'C:\Users\Jennifer\Documents\MATLAB\Stimuli\DRCs'
filename = 'DRC002';

% load filter
filtName = 'C:\Users\Jennifer\Documents\filters\180816_blueAcuteSpeaker_NIDAQ_3k-80k_fs200k.mat';
load(filtName);

% parameters
params.seed = 89; %13
params.fs = 200000; %sampling rate
params.nTones = 50; %number of tones
params.freqs = logspace(log10(5000),log10(40e3),params.nTones);
params.mu = 50; %Average sound level
params.sd = 20;
params.rampDuration = .001; % in seconds
params.chordDuration = .02; % in seconds (total chord duration + ramp duration)
params.blockDuration = 300; % in seconds
params.baseAmplitude = .1; % Correction for calibration
params.nBlocks = 1;
params.totalDuration = params.nBlocks * params.blockDuration;
params.nFiles = 8;
params.filter = filtName;
fs = params.fs;

laserON = 0.5;
laserdur = 0.25;
evpulse = 0.005;

ev1 = zeros(fs,1);
ev1(1:evpulse*fs) = 0.5;
laser1 = zeros(fs,1);
laser1(laserON*fs:(laserON+laserdur)*fs) = 0.5;

ev = repmat(ev1,params.blockDuration,1);
laser = repmat(laser1,params.blockDuration,1);

% make some noise
for n = 1:params.nFiles
    fprintf('file %02d... ', n)
    fprintf('making waveform... ')
    [stim amps dbs] = makeContrastBlocks(params,params.nBlocks, 'uniform','spline'); 
    
    %filter
    fprintf('filtering waveform... ')
    stimFilt = conv(stim,FILT,'same');
    STIM = [stimFilt' ev laser];
    
    fprintf('saving waveform... \n')
    save([filename '_180816_blueAcuteSpeaker' '-' sprintf('%02d',n)],'params','amps','dbs');
    %audiowrite([filename '-' sprintf('%02d',n) '.wav'],stim, params.fs);
    audiowrite([filename '_180816_blueAcuteSpeaker' '-' sprintf('%02d',n) '.wav'],STIM, params.fs);
end

%soundsc([stim(1:2*fs)], fs)
%spectrogram([stim(1:2*fs)],256,64,params.freqs,fs,'yaxis');