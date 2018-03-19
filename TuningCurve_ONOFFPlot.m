function TuningCurve_ONOFFPlot(TCon, TCoff)
% A not very elegant way to plot both the laser on and laser off tuning
% curves on a single plot.
%
% Inputs:
%   TCon: Tuning curve output from TuningCurve_Compute.m running laser ON
%         trials
%   TCon: Tuning curve output from TuningCurve_Compute.m running laser OFF
%         trials
%   TITLE: String, title of figure




subplot(2, 4, 1);
OneChannelPlot(TCon.FRmat', TCon.t, TCon.F', 3, 'Time (s)', 'Frequency (Hz)', 'Firing rate');
if length(TCon.CellInfo) == 4
    TCon.CellInfo = [TCon.CellInfo 0 8];
end
title('Laser On');

subplot(2, 4, 3);
if (length(TCon.A) == 1)
    plot(TCon.F, TCon.TCmat{1});
    box off, axis tight;
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC 0-.08 s');
    box off;
else
    OneChannelPlot(TCon.TCmat{1}, TCon.F', TCon.A', 1,  'Frequency (Hz)','Amplitude', 'TC 0-.08 s');
    colorbar;
end

subplot(2, 4, 5);
if (length(TCon.A) == 1)
    plot(TCon.F, TCon.TCmat{2});
    box off, axis tight;
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC .1-.15 s')
    box off;
else
    OneChannelPlot(TCon.TCmat{2}, TCon.F', TCon.A', 1,  'Frequency (Hz)','Amplitude', 'TC .1-.15 s');
end

subplot(2, 4, 7);
if (length(TCon.A) == 1)
    plot(TCon.F, TCon.TCmat{3});
    box off;
    axis tight
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC .2-.25 s')
    box off;
else
    OneChannelPlot(TCon.TCmat{3},  TCon.F', TCon.A', 1,  'Frequency (Hz)','Amplitude', 'TC .2-.25 s');
end

%Plot laser OFF data
subplot(2, 4, 2);
OneChannelPlot(TCoff.FRmat', TCoff.t, TCoff.F', 3, 'Time (s)', 'Frequency (Hz)', 'Firing rate');
title('Laser Off');

subplot(2, 4, 4);
if (length(TCoff.A) == 1)
    plot(TCoff.F, TCoff.TCmat{1});
    box off, axis tight;
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC 0-.08 s');
    box off;
else
    OneChannelPlot(TCoff.TCmat{1}, TCoff.F', TCoff.A', 1,  'Frequency (Hz)','Amplitude', 'TC 0-.08 s');
    colorbar;
end

subplot(2, 4, 6);
if (length(TCoff.A) == 1)
    plot(TCoff.F, TCoff.TCmat{2});
    box off, axis tight;
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC .1-.15 s')
    box off;
else
    OneChannelPlot(TCoff.TCmat{2}, TCoff.F', TCoff.A', 1,  'Frequency (Hz)','Amplitude', 'TC .1-.15 s');
end

subplot(2, 4, 8);
if (length(TCoff.A) == 1)
    plot(TCoff.F, TCoff.TCmat{3});
    box off;
    axis tight
    set(gca, 'xscale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    title('TC .2-.25 s')
    box off;
else
    OneChannelPlot(TCoff.TCmat{3},  TCoff.F', TCoff.A', 1,  'Frequency (Hz)','Amplitude', 'TC .2-.25 s');
end

%Set lims for both to be the same
fmax1 = max([max(TCoff.TC_FR(:)) max(TCon.TC_FR(:))]);
fmax2 = max([max(TCoff.TCmat{1}(:)) max(TCon.TCmat{1}(:))]);
for n = 1:2
    subplot(2, 4, n);
    if fmax1 > 0;
        set(gca, 'cLim', [0 fmax1]);
    end
end
for n = 3:8
    subplot(2, 4, n);
    if fmax2 > 0;
        if length(TCoff.A) > 1
            set(gca, 'cLim', [0 fmax2]);
        else
            set(gca,'YLim',[0 fmax2]);
        end
    end
end
suptitle(num2str(TCon.CellInfo));