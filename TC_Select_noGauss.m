function [LOCSon,LOCSoff] = TC_Select_noGauss(MNum, h1, nfreq, amps, fig)
% This function is meant to locate nfreq frequencies with highest firing rates
% for each cell
%INPUTS:
%   MNum = mouse number
%   h1 = List of files with extra 0 ammended, but not including filter
%   number
%   nfreq = number of frequencies of interest
%   amps = amplitude indices you want to average across (array)
%   fig = 1 to display figures, 0 to not display figures.
%
%OUTPUTS:
%   LOCSon: indices of top frequency locations for laser on condition
%   LOCSoff: indices of top frequency locations for laser off condition

LOCSon = NaN(size(h1,1),nfreq);
LOCSoff = NaN(size(h1,1),nfreq);
nAmps = length(amps);

cd(['D:\Spikes\M' num2str(MNum) '\TCs']);
for v = 1:size(h1, 1); 
    
    %Load and smooth tuning curves:
    load(['data\' h1(v,:)]);
    
    if exist('TC1','var') %Added this statement because partway through this project I started renaming this variable.
        for u = 1:nAmps
            T1(:, u) = SmoothGaus(TC1.TCmat{1}(:, amps(u)), 3);
        end

        for u = 1:nAmps
            T2(:, u) = SmoothGaus(TC2.TCmat{1}(:, amps(u)), 3);
        end
    else
        for u = 1:nAmps
            T1(:, u) = SmoothGaus(TCon.TCmat{1}(:, amps(u)), 3);
        end

        for u = 1:nAmps
            T2(:, u) = SmoothGaus(TCoff.TCmat{1}(:, amps(u)), 3);
        end
    end
    
    %Take mean across top desired amplitudes
    laser = mean(T1(:,1:nAmps),2);
    nolaser = mean(T2(:,1:nAmps),2);
    
    %Find indices for top nfreq frequencies in nolaser condition
    [~,IDXoff] = sort(nolaser,'descend');
    [~,IDXon] = sort(laser,'descend');
    
    LOCSoff(v,:) = IDXoff(1:nfreq)';
    LOCSon(v,:) = IDXon(1:nfreq)';
    
   
   %Figures dispalying fits and comparing laser on vs laser off tuning
   if fig == 1
        figure; plot(TC2.F, laser, 'r');
        hold on; plot(TC2.F, nolaser, 'k');
        set(gca, 'xscale', 'log');        
        box off; hold off;
        title(h1(v,9:5+q-1))
   end
   
    
end


