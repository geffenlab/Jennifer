function [GAUSSparams, LOCS] = TC_Select5(MNum, h1, fig)
% This function is meant to  help select for good cells based on having
% single peaked tuning curves. It also finds best frequency and tuning
% curve width.
%INPUTS:
%   MNum = mouse number
%   h1 = List of files with extra 0 ammended, but not including filter
%   number
%   filt = laser filter value (which will be present as last number in
%   file name, but not in h1)
%   fig = 1 to display figures, 0 to not display figures.

%OUTPUTS:
%   GAUSSparams = parameters of gaussian fit, including R^2 value
%   LOCS: index of best frequency location

cd(['D:\Spikes\M' num2str(MNum) '\TCs']);

for v = 1:size(h1, 1); 
    
    %Load and smooth tuning curves:
    load(['data\' h1(v,:)]);
    
    for u = 1:8
        T1(:, u) = SmoothGaus(TC1.TCmat{1}(:, u), 3);
    end
    
    for u = 1:8
        T2(:, u) = SmoothGaus(TC2.TCmat{1}(:, u), 3);
    end
    
    %Take mean across top 3 amplitudes
    laser = mean(T1(:,6:8),2);
    nolaser = mean(T2(:,6:8),2);


%Calculate gaussian fit:
    x = TC1.F';
    norm_func = @(b,x) b(1) .* exp(-((x - b(2))/b(3)).^2) + b(4);
    initguess = [5 15000 5000 0]; %Inital guess of paramters [amplitude, center freq, sdev, baseline shift]
    coeffval_LaserOFF = nlinfit(x,nolaser,norm_func,initguess); 
    coeffval_LaserON = nlinfit(x,laser,norm_func,initguess);
    y_OFF = norm_func(coeffval_LaserOFF,x); %Evaluate function with found parameters
    y_ON = norm_func(coeffval_LaserON,x); %Evaluate function with found parameters

    %Start pulling out paramters:
    rsq_ON = Rsquared(laser,y_ON); %Calculate R^2 value 
    amp_ON = coeffval_LaserON(1); %Max amplitude
    bf_ON = coeffval_LaserON(2); %Best frequency
    f1_ON = coeffval_LaserON(2)-coeffval_LaserON(3)*sqrt(log(2)); %lower frequency at 1/2 max amplitude
    f2_ON = coeffval_LaserON(2)+coeffval_LaserON(3)*sqrt(log(2)); %higher frequensy at 1/2 max amplitude
    bw_ON = f2_ON-f1_ON; % 1/2 max bandwidth
    bw_octave_ON = log2(f2_ON/f1_ON);
    
    rsq_OFF = Rsquared(nolaser,y_OFF);
    amp_OFF = coeffval_LaserOFF(1);
    bf_OFF = coeffval_LaserOFF(2);
    f1_OFF = coeffval_LaserOFF(2)-coeffval_LaserOFF(3)*sqrt(log(2)); %lower frequency at 1/2 max amplitude
    f2_OFF = coeffval_LaserOFF(2)+coeffval_LaserOFF(3)*sqrt(log(2)); %higher frequency at 1/2 max amplitude
    bw_OFF = f2_OFF-f1_OFF; % 1/2 max bandwidth
    bw_octave_OFF = log2(f2_OFF/f1_OFF);

    GAUSSparams.rsqON(v) = rsq_ON;
    GAUSSparams.rsqOFF(v) = rsq_OFF;
    GAUSSparams.BestFreqON(v) = bf_ON;
    GAUSSparams.BestFreqOFF(v) = bf_OFF;
    GAUSSparams.WidthONoct(v) = bw_octave_ON;
    GAUSSparams.WidthOFFoct(v) = bw_octave_OFF;
    GAUSSparams.WidthON(v) = bw_ON;
    GAUSSparams.WidthOFF(v) = bw_OFF;
    
   
   %Figures dispalying fits and comparing laser on vs laser off tuning
   if fig == 1
        figure; subplot(2,2,[1 2]); plot(TC2.F, laser, 'r');
        hold on; plot(TC2.F, nolaser, 'k');
        set(gca, 'xscale', 'log');        
        box off; hold off;
        title([h1(v,9:5+q-1) '-' num2str(filt)])

        subplot(2,2,3); plot(x,y_OFF,'b')
        hold on; plot(x,nolaser,'k');
        set(gca, 'xscale', 'log'); 
        box off; hold off;
        title(['Laser OFF fit; R^2 = ' num2str(rsq_OFF)])
        
        subplot(2,2,4);plot(x,y_ON,'b')
        hold on; plot(x,laser,'r')
        set(gca, 'xscale', 'log'); 
        box off; hold off
        title(['Laser ON fit; R^2 = ' num2str(rsq_ON)])
        
        %print('-djpeg','-r500',[h1(v,9:5+q-1) '-' num2str(filt)]);
   end
   
   %Index of best frequency
   if rsq_ON > 0.65 || rsq_OFF > 0.65
       [~,LOCS.laser(v)] = min(abs(x - bf_ON));
       [~,LOCS.nolaser(v)] = min(abs(x - bf_OFF));
   else LOCS.laser(v) = NaN; LOCS.nolaser(v) = NaN;
   end
    
end

%SOMETIMES makes big mistakes with curves, NEED AN ALTERNATIVE FOR FINDING BF


