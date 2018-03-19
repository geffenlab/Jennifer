function TuningCurve_Plot(TC, TITLE)
%Plot tuning curve
%
% Inputs:
%   TC: Tuning curve output from TuningCurve_Compute.m
%   TITLE: String, title of figure


subplot(2, 2, 1);
suptitle(TITLE);
OneChannelPlot(TC.FRmat', TC.t, TC.F', 3, 'Time (s)', 'Frequency (Hz)', 'Firing rate');
if length(TC.CellInfo) == 4
    TC.CellInfo = [TC.CellInfo 0 8];
end
title(num2str(TC.CellInfo));

subplot(2, 2, 2);
if (length(TC.A) == 1)
    plot(TC.F, TC.TCmat{1});
    box off, axis tight;
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC 0-.08 s');
    box off;
else
    OneChannelPlot(TC.TCmat{1}, TC.F', TC.A', 1,  'Frequency (Hz)','Amplitude', 'TC 0-.08 s');
    colorbar;
end

subplot(2, 2, 3);
if (length(TC.A) == 1)
    plot(TC.F, TC.TCmat{2});
    box off, axis tight;
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC .1-.15 s')
    box off;
else
    OneChannelPlot(TC.TCmat{2}, TC.F, TC.A, 1,  'Frequency (Hz)','Amplitude', 'TC .1-.15 s');
end

subplot(2, 2, 4);
if (length(TC.A) == 1)
    plot(TC.F, TC.TCmat{3});
    box off;
    axis tight
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC .2-.25 s')
    box off;
else
    OneChannelPlot(TC.TCmat{3},  TC.F, TC.A, 1,  'Frequency (Hz)','Amplitude', 'TC .2-.25 s');
end
