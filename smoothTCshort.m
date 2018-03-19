function [FR_NOLASER, FR_LASER, CellQ] = smoothTCshort(MNum,file,nfreq,amps)
%Function to smooth tuning curve firing rate across top nfreq and
%designated amplitudes.
%
% Parameters
% ----------
%   MNum: integer
%       Mouse identifying number
%   file: string
%       Contains file name string of format TC001_LEFT-0X-MXXXX.....mat
%   nfreq: integer vector
%       Number of frequencies to use in smoothed time course
%   amps: integer vector
%       Number of amplitudes to use in smoothed time course
%       
%
% Returns
% -------
%   FR_NOLASER: vector
%       Contains smoothed firing time course for tuning curve stim for no laser condition
%   FR_LASER: vector
%       Contains smoothed firing time course for tuning curve stim for laser condition 
%   CellQ: integer
%       Cell quality


    %Pre-allocate variables
    allCELL = [];
    FR_LASER = []; 
    FR_NOLASER = [];
    %LOCSoff = [];
    %LOCSon = [];
    CellQ = [];
    all_LASERSPIKE2 = []; 
    all_NOLASERSPIKE2 = [];
    LASER_SPIKE = [];
    NOLASER_SPIKE = [];
    SpikeDataLaser = [];        
    SpikeDataNoLaser = [];
    LASER_SPIKE2 = {};
    NOLASER_SPIKE2 = {};
    
    
    %Step 1: Must find the frequencies to use.   
    %nfreq = 7; %Number of frequencies to use in smoothes response
    [LOCSon,LOCSoff] = TC_Select_noGauss(MNum,file,nfreq,0); %Need to adjust this code for new stimulus
    
    %Step 2: Separate spikes by tone trial (separated by trial, frequency,
    %amplitude, and laser vs no laser)
    afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
    Win = 0.24; %Window (in seconds) around start of tone to look for spikes

    %Separate out spikes based on if during laser trials vs no laser trials
    load(['D:\Spikes\M' num2str(MNum) '\SpikeMat\'  file '.mat']); 
    CellQ = CellInfo(6);

    [SpkTime_Laser, SpkTime_NoLaser, ~, ~] = SpikeTime(afo5,SpikeData,nRep,Win);

    %Select those in bins around top frequencies and designated
    %amplitudes
    LASER_SPIKE = SpkTime_Laser(LOCSon,amps,:);
    NOLASER_SPIKE = SpkTime_NoLaser(LOCSoff,amps,:);

    %Step 3: Calculate and smooth firing rate for laser and no laser conditions
    for w = 1:numel(LASER_SPIKE)
        LASER_SPIKE2{w} = [sort(LASER_SPIKE{w}); w*ones(1,length(LASER_SPIKE{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
        SpikeDataLaser = [SpikeDataLaser LASER_SPIKE2{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
    end
    FR_LASER = smoothFRx4(SpikeDataLaser,numel(LASER_SPIKE),0.001,[-Win Win],5);
    all_LASERSPIKE2 = LASER_SPIKE2;

    for w = 1:numel(NOLASER_SPIKE)
        NOLASER_SPIKE2{w} = [sort(NOLASER_SPIKE{w}); w*ones(1,length(NOLASER_SPIKE{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
        SpikeDataNoLaser = [SpikeDataNoLaser NOLASER_SPIKE2{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
    end
    FR_NOLASER = smoothFRx4(SpikeDataNoLaser,numel(NOLASER_SPIKE),0.001,[-Win, Win],5);
    all_NOLASERSPIKE2 = NOLASER_SPIKE2;

        
