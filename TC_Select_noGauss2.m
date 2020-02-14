function [LOCSon,LOCSoff] = TC_Select_noGauss2(filename, freqidx, amps,n)
% This function is meant to locate nfreq frequencies with highest firing rates
% for each cell
%INPUTS:
%   MNum = mouse number
%   h1 = File name containing TC1/TC2 or TCon/TCoff files
%   freqidx = frequency indexes you want to use (1 being most preferred)
%   amps = amplitude indices you want to average across (array)
%   fig = 1 to display figures, 0 to not display figures.
%
%OUTPUTS:
%   LOCSon: indices of top frequency locations for laser on condition
%   LOCSoff: indices of top frequency locations for laser off condition

nfreq = length(freqidx);
nAmp = length(amps);

LOCSon = NaN(1,nfreq);
LOCSoff = NaN(1,nfreq);

cd(['D:\Spikes\' filename(15:19) '\TCs\data']);
load(['TC003-' filename(15:end-4) '_laser-0' num2str(n) '.mat'])
%for v = 1:size(h1, 1); 
    
    %Load and smooth tuning curves:
    %load(['data\' h1(v,:)]);
    
    if exist('TC1','var') %Added this statement because partway through this project I started renaming this variable.
        for u = 1:nAmp
            T1(:, u) = SmoothGaus(TC1.TCmat{1}(:, amps(u)), 3);
        end

        for u = 1:nAmp
            T2(:, u) = SmoothGaus(TC2.TCmat{1}(:, amps(u)), 3);
        end
    else
        for u = 1:nAmp
            T1(:, u) = SmoothGaus(TCon.TCmat{1}(:, amps(u)), 3);
        end

        for u = 1:nAmp
            T2(:, u) = SmoothGaus(TCoff.TCmat{1}(:, amps(u)), 3);
        end
    end
    
    %Take mean across top desired amplitudes
    laser = mean(T1(:,1:nAmp),2);
    nolaser = mean(T2(:,1:nAmp),2);
    
    %Find indices for top nfreq frequencies in nolaser condition
    [~,IDXoff] = sort(nolaser,'descend');
    [~,IDXon] = sort(laser,'descend');
    
    LOCSoff = IDXoff(freqidx)';
    LOCSon = IDXon(freqidx)';
    
   
    
end


