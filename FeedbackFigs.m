FigOut = 'C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\Manuscripts\Feedback\Figures';
tunedata = 'C:\Users\Jennifer\Documents\MATLAB\TuningCurveAnalysis\Feedback ChR2 in IC.mat';
load(tunedata);

GOODCELLall = vertcat(GOODCELL{:});


fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k','defaultAxesFontSize',16);


%% ************************************************************************
%  *****                       1. SILENCE + LASER                     *****
%  ************************************************************************


%%%%%%%%%%%%% I. EXAMPLE TRACE %%%%%%%%%%%%%
% Plot Raster and smoothed FR for both 25 ms and 250 ms stimuli (these were
% the two that significantly affect FR)

EX = {'M3222-2-9-1-1'; 'M3222-3-8-2-1'}; %Example units
for i = 1:length(EX)
    load(['R452-' EX{i} '.mat']); %25 ms stim
    smoothedFR = smoothFRx4(SpikeData(3,:),nRep,0.001,[0 1],5);
    
    subplot(2,2,1); %Raster
    line([SpikeData(3,:); SpikeData(3,:)], [SpikeData(4,:); SpikeData(4,:)+1],'color','k','linewidth',1.5)
    hold on; line([0.1 0.1],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--');
         line([0.125 0.125],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--'); 
    set(gca,'TickDir','out','XTick',[0:0.2:1],'XTickLabel',[0:0.2:1])
    axis tight; box off;
    ylabel('Trial #');
    
    subplot(2,2,3); %Smooth FR
    plot(smoothedFR,'k','linewidth',2)
    set(gca,'TickDir','out','XTick',[0:200:1000],'XTickLabel',[0:0.2:1])
    axis tight; box off;
    xlabel('Time (s)'); ylabel('Firing Rate (Hz)')
    
    load(['R453-' EX{i} '.mat']); %250 ms stim
    smoothedFR = smoothFRx4(SpikeData(3,:),nRep,0.001,[0 1],5);

    subplot(2,2,2); %Raster
    line([SpikeData(3,:); SpikeData(3,:)], [SpikeData(4,:); SpikeData(4,:)+1],'color','k','linewidth',1.5)
    hold on; line([0.1 0.1],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--');
         line([0.35 0.35],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--'); 
    set(gca,'TickDir','out','XTick',[0:0.2:1],'XTickLabel',[0:0.2:1])
    axis tight; box off;
    ylabel('Trial #');
    
    subplot(2,2,4); %Smooth FR
    plot(smoothedFR,'k','linewidth',2)
    set(gca,'TickDir','out','XTick',[0:200:1000],'XTickLabel',[0:0.2:1])
    axis tight; box off;
    xlabel('Time (s)'); ylabel('Firing Rate (Hz)')
    suptitle(EX{i});
    
    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperUnits','points');
    set(gcf,'PaperSize',[1600 800]);
    set(gcf,'Position',[0 0 1600 800]);

end

%%%%%%%%%%%%% II. PLOT LASER ON VS LASER OFF RESPONSE %%%%%%%%%%%%%

Stim = {'R450';'R451';'R452';'R453'};
onset = 100; %Laser onset (ms)
dur = [1 5 25 250]; %Laser duration (ms) 
SilenceON = NaN(size(Qcell,1),length(Stim));
SilenceOFF = NaN(size(Qcell,1),length(Stim));

% Use all cells (multiunits + single) to look for laser affect
temp = vertcat(allCELL{:});
Qcell = temp(find(horzcat(CellQ{:}) > 0 & horzcat(CellQ{:}) < 5),:);

for u = 1:size(Qcell,1)
    
    for i = 1:length(Stim)
        q = find(Qcell(u,:) == '.');
        if exist([Stim{i} '-' Qcell(u,7:q - 10) '.mat'])
            load([Stim{i} '-' Qcell(u,7:q - 10) '.mat']);
        end
        if ~isempty(SpikeData)
            smoothedFR = smoothFRx4(SpikeData(3,:),nRep,0.001,[0 1],5);
            FRmeanON = mean(smoothedFR(onset:onset+dur(i)));
            FRmeanOFF = mean(smoothedFR(1:onset-1));
            if FRmeanOFF > 0.5 % Low spontaneous firing rates can lead to outliers
                SilenceON(u,i) = FRmeanON;
                SilenceOFF(u,i) = FRmeanOFF;
            end            
        end
    end
end

%Plot laser ON vs laser OFF and then mean for all 4 laser durations
colour = [0.2 0.2 0.2; 0.4 0.4 0.4; 0.6 0.6 0.6; 0.8 0.8 0.8];
subplot(4,2,[1 2 3 4]);
for i = 1:length(Stim)
    hold on;
    scatter(SilenceOFF(:,i), SilenceON(:,i),25,'filled','MarkerFaceColor', colour(i,:))
    [p_silence(i),~] = signrank(SilenceOFF(:,i), SilenceON(:,i));
end
max1 = nanmax(SilenceOFF(:));
max2 = nanmax(SilenceON(:));
line([0 nanmax(max1,max2)],[0 nanmax(max1,max2)],'linestyle','--','color','k');
legend(['1 ms, p = ' num2str(p_silence(1))] ,['5 ms, p = ' num2str(p_silence(2))],['25 ms, p = ' num2str(p_silence(3))],['250 ms, p = ' num2str(p_silence(4))],'Location','best')
axis square; axis tight
xlabel('Firing Rate OFF (Hz)'); ylabel('Firing Rate ON (Hz)');
set(gca,'TickDir','out')

for i = 1:length(Stim)
    subplot(4,2,i+4);
    bar([nanmean(SilenceOFF(:,i)) nanmean(SilenceON(:,i))],0.5,'EdgeColor','none','FaceColor',colour(i,:));
     hold on; errorbar([nanmean(SilenceOFF(:,i)) nanmean(SilenceON(:,i))],[nanstd(SilenceOFF(:,i))...
        ./sqrt(sum(isnan(SilenceOFF(:,i)))) nanstd(SilenceON(:,i))./sqrt(sum(isnan(SilenceON(:,i))))],'k','LineStyle','none','LineWidth',2)
    box off; ylabel('Firing Rate (Hz)')
    set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square

end
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 1600]);
set(gcf,'Position',[0 0 800 1600]);

%% ************************************************************************
%  *****                       2. CLICK FIGURES                       *****
%  ************************************************************************

ClickOFF = [0 2 4]; %Onset of clicks in laser off condition
ClickON = [1 3 5]; %Onset of clicks in laser on condition
ClickDUR = 0.05; %Window (in seconds) for single click duration
LaserON = [0.75 2.75 4.75]; %Laser onset (seconds)
LaserDUR = 1; %Laser duration (seconds)

%%%%%%%%%%%%% I. EXAMPLE TRACE %%%%%%%%%%%%%
EX = {'M3222-2-8-1-1'};
load(['R400-' EX{1} '.mat']);

subplot(2,1,1);
line([SpikeData(3,:); SpikeData(3,:)], [SpikeData(4,:); SpikeData(4,:)+1],'color','k','linewidth',1.5)
hold on; line([LaserON(1) LaserON(1)],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--');
line([LaserON(1)+LaserDUR LaserON(1)+LaserDUR],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--'); 
line([LaserON(2) LaserON(2)],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--');
line([LaserON(2)+LaserDUR LaserON(2)+LaserDUR],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--'); 
line([LaserON(3) LaserON(3)],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--');
line([LaserON(3)+LaserDUR LaserON(3)+LaserDUR],[0 nRep+1],'color','m','linewidth',1.5,'linestyle','--'); 
box off; axis tight;
set(gca,'TickDir','out','XTick',[0:1:6],'XTickLabel',[0:1:6]);
ylabel('Trial #')

smoothedFR = smoothFRx4(SpikeData(3,:),nRep,0.001,[0 6],5);
subplot(2,1,2); plot(smoothedFR,'k','linewidth',2);
box off; axis tight; 
set(gca,'TickDir','out','XTick',[0:1000:6000],'XTickLabel',[0:1:6]);
xlabel('Time (s)'); ylabel('Firing Rate (Hz)');
suptitle(EX{1});

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 800]);
set(gcf,'Position',[0 0 1600 800]);



%%%%%%%%%%%%% II. Laser ON vs Laser OFF RESPONSE %%%%%%%%%%%%%

% Use all cells (multiunits + single) to look for laser affect
temp = vertcat(allCELL{:});
Qcell = temp(horzcat(CellQ{:}) > 0 & horzcat(CellQ{:}) < 5,:);


%Pre-allocate variables
spontON = NaN(1,size(Qcell,1));
spontOFF = NaN(1,size(Qcell,1));
ClickOFF = NaN(1,size(Qcell,1));
ClickON = NaN(1,size(Qcell,1));

for u = 1:size(Qcell,1)
    q = find(Qcell(u,:) == '.');
    if exist(['R400-' Qcell(u,7:q - 10) '.mat'],'file')
        load(['R400-' Qcell(u,7:q - 10) '.mat']);
    end
    if ~isempty(SpikeData)
        idxOFF1 = find(SpikeData(3,:) > 0 & SpikeData(3,:) < 1);
        idxOFF2 = find(SpikeData(3,:) > 2 & SpikeData(3,:) < 3);
        idxOFF3 = find(SpikeData(3,:) > 4 & SpikeData(3,:) < 5);
        SpikesOFF = SpikeData(:,[idxOFF1 idxOFF2 idxOFF3]);
        smoothOFF = smoothFRx4(mod(SpikesOFF(3,:),1),nRep*3,0.001,[0 1],5);
        if mean(smoothOFF(1:50)) > 2*std(smoothOFF(650:700)) %Shows click response
            spontON(u) = mean(smoothOFF(750:end));
            ClickOFF(u) = mean(smoothOFF([1:50 100:150 200:250 300:350 400:450 500:550]));

            idxON1 = find(SpikeData(3,:) > 1 & SpikeData(3,:) < 2);
            idxON2 = find(SpikeData(3,:) > 3 & SpikeData(3,:) < 4);
            idxON3 = find(SpikeData(3,:) > 5 & SpikeData(3,:) < 6);
            SpikesON = SpikeData(:,[idxON1 idxON2 idxON3]);
            smoothON = smoothFRx4(mod(SpikesON(3,:),1),nRep*3,0.001,[0 1],5);
            spontOFF(u) = mean(smoothON(750:end));
            ClickON(u) = mean(smoothON([1:50 100:150 200:250 300:350 400:450 500:550]));
        end
        
    end
    
end
[p_click, h_click] = signrank(ClickOFF, ClickON);
figure; subplot(1,2,1);
scatter(ClickOFF,ClickON,25,'filled','MarkerFaceColor',[0.2 0.2 0.2])
max1 = nanmax(ClickOFF(:));
max2 = nanmax(ClickON(:));
line([0 nanmax(max1,max2)],[0 nanmax(max1,max2)],'linestyle','--','color','k');
axis square; axis tight
set(gca,'TickDir','out');
xlabel('Firing Rate OFF (Hz)'); ylabel('Firing Rate ON (Hz)'); title(['p = ' num2str(p_click)]);
subplot(1,2,2);
 bar([nanmean(ClickOFF) nanmean(ClickON)],0.5,'EdgeColor','none','FaceColor',[0.2 0.2 0.2]);
     hold on; errorbar([nanmean(ClickOFF) nanmean(ClickON)],[nanstd(ClickOFF)...
        ./sqrt(sum(isnan(ClickOFF))) nanstd(ClickON)./sqrt(sum(isnan(ClickON)))],'k','LineStyle','none','LineWidth',2)
box off; ylabel('Firing Rate (Hz)'); title(['N = ' num2str(sum(isnan(ClickOFF)))]);
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
suptitle('Click Response')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1000 600]);
set(gcf,'Position',[0 0 1000 600]);


%%%%%%%%%%%%% III. CHANGE IN SPONT VS CHANGE IN CLICK RESPONSE %%%%%%%%%%%%%
%data calculated in previous plot section
max1 = nanmax([spontON - spontOFF]);
max2 = nanmax(ClickON - ClickOFF);
min1 = nanmin(spontON - spontOFF);
min2 = nanmin(ClickON - ClickOFF);

figure;
scatter([spontON - spontOFF],[ClickON - ClickOFF],25,'filled')
hold on; line([min(min1, min2) max(max1, max2)], [min(min1, min2) max(max1, max2)],'linestyle','--','color','k')
axis tight; box off; axis square
set(gca,'TickDir','out')
xlabel('delta Spontaneous activity (ON - OFF)'); ylabel('delta Click response (ON - OFF)')



%% ************************************************************************
%  *****                       3. TUNING FIGURES                      *****
%  ************************************************************************

%%%%%%%%%%%%% I. EXAMPLE TCs AND TIMECOURSES %%%%%%%%%%%%%

tLaserON = [0.4 0.48 0.508];
tLaserOFF = [0.65 0.73 0.758];
EX = {'M3222-1-5-3-1';'M3222-1-7-3-1'; 'M3223-1-8-3-1'; 'M3226-3-6-1-1'};
for ex = 1:length(EX)
    for i = 1:length(tLaserON)
        cd('D:\Spikes')
        load(['TC003-' EX{ex} '_laser-0' num2str(i)]);
        subplot(2, 3, i);
        tt1 = find(TCon.t>=.2, 1, 'first');
        tt5 = find(TCon.t>=.8, 1, 'first');
        fr1 = SmoothGaus(mean(TCon.TC_FR, 1), 3);
        fr2 = SmoothGaus(mean(TCoff.TC_FR, 1), 3);
        plot(TCoff.t(tt1:tt5), fr1(tt1:tt5), 'k','linewidth',2);
        hold on;
        plot(TCoff.t(tt1:tt5), fr2(tt1:tt5), 'r','linewidth',2);
        ylabel('Firing rate (Hz)'); xlabel('Time (s)');
        line([tLaserON(i) tLaserOFF(i); tLaserON(i) tLaserOFF(i)], [0 0; max(fr1) max(fr1)], 'color', 'm','linestyle','--');
        line([0.5 0.55; 0.5 0.55], [0 0; max(fr1) max(fr1)], 'color', 'k','linestyle','--');
        box off; axis tight;
        hold off;
        set(gca,'TickDir','out')

        Ton = SmoothGaus(TCon.TCmat{1}, 3);
        Toff = SmoothGaus(TCoff.TCmat{1}, 3);


        subplot(2,3,3+i)
        plot(TCoff.F, Ton, 'r','linewidth',2);
        hold on; plot(TCoff.F, Toff, 'k', 'linewidth',2);
        set(gca, 'xscale', 'log'); 
        xlabel('Frequency (Hz)')
        ylabel('Firing rate (Hz)')
        box off; axis tight;
        hold off;
        set(gca,'TickDir','out')
    end
suptitle(EX{ex})  
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[2000 1000]);
set(gcf,'Position',[0 0 2000 1000]);

cd(FigOut)
print(['Example' num2str(ex)],'-dpdf','-r400')
end

%%
%%%%%%%%%%%%% II. SPARSENESS %%%%%%%%%%%%%

%Load tuning curve parameters
load('D:\Code\TuningCurve\TC003_170815LEFT_stimInfo');
fList = stimInfo.index(:,1);
aList = stimInfo.attenuations;
StimOrder = [fList(stimInfo.order)'; aList(ones(length(stimInfo.order),1))'];
trialdur = (stimInfo.tDur+stimInfo.ITI)/1000;
dur = trialdur*length(stimInfo.order);
Time = 0:trialdur:dur; %Each repetition is 400 seconds long, each trial is 500ms long
Time = Time(1:end-1);

StimOrder_Laser = [Time(2:2:end); StimOrder(:,2:2:end)];
StimOrder_NoLaser = [Time(1:2:end); StimOrder(:,1:2:end)];

nFreq = 50;
amps = 1;

TITLES = {'Response Mag DOWN';'Response Mag UP';'spont DOWN';'spont UP'};


cellidx{1} = GOODCELLall(SigCellmagDOWN,:);
cellidx{2} = GOODCELLall(SigCellmagUP,:);
cellidx{3} = GOODCELLall(SigCellspontDOWN,:);
cellidx{4} = GOODCELLall(SigCellspontUP,:);

for ii = 1%:4
    h1 = cellidx{ii};
    numcell = size(h1,1);
    SidxON = nan(3,numcell);
    SidxOFF = nan(3,numcell);
    FR_LASER_ALL = [];
    FR_NOLASER_ALL = [];

    for n = 1:3


        Win = 0.24; %Window (in seconds) around start of tone to look for spikes

        SpkTime_LaserAll = cell(size(h1,1),nFreq);
        SpkTime_NoLaserAll = cell(size(h1,1),nFreq);

        for v = 1:size(h1,1) 

            match = strfind(h1(v,:), '_');
            q = match(end);
            load(['D:\Spikes\M' h1(v,8:11) '\SpikeMat\TC003_LEFT-0' num2str(n) h1(v,6:q-1) '.mat']); 

            SpkTime_NoLaser = SpikeTime(StimOrder_NoLaser,SpikeData,nRep, Win);
            SpkTime_Laser = SpikeTime(StimOrder_Laser,SpikeData,nRep, Win);

            for j = 1:nFreq
                SpkTime_LaserAll{v,j} = SpkTime_Laser(j,amps,:);
                SpkTime_NoLaserAll{v,j} = SpkTime_NoLaser(j,amps,:);
            end
        end


        %Calculate and smooth firing rate (NOTE: difference from following
        %section is that this includes all frequencies and all amplitudes)
        SpikeDataLaserAll = cell(size(h1,1),nFreq);
        SpikeDataNoLaserAll = cell(size(h1,1),nFreq);
        clear LASER_SPIKE_ALL NOLASER_SPIKE_ALL
        %For laser trials
        for j = 1:nFreq
            for v = 1:size(h1,1)
                for w = 1:numel(SpkTime_LaserAll{v})
                    LASER_SPIKE_ALL{v,j}{w} = [sort(SpkTime_LaserAll{v,j}{w}); w*ones(1,length(SpkTime_LaserAll{v,j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                    SpikeDataLaserAll{v,j} = [SpikeDataLaserAll{v,j} LASER_SPIKE_ALL{v,j}{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
                end
                FR_LASER_ALL(j,v,:) = smoothFRx4(SpikeDataLaserAll{v,j},numel(LASER_SPIKE_ALL{v,j})*stimInfo.repeats,0.001,[-Win Win],5);

                %For no laser trials
                for w = 1:numel(SpkTime_NoLaserAll{v})
                    NOLASER_SPIKE_ALL{v,j}{w} = [sort(SpkTime_NoLaserAll{v,j}{w}); w*ones(1,length(SpkTime_NoLaserAll{v,j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                    SpikeDataNoLaserAll{v,j} = [SpikeDataNoLaserAll{v,j} NOLASER_SPIKE_ALL{v,j}{w}];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
                end
                FR_NOLASER_ALL(j,v,:) = smoothFRx4(SpikeDataNoLaserAll{v,j},numel(NOLASER_SPIKE_ALL{v,j})*stimInfo.repeats,0.001,[-Win Win],5);

            end
        end


        %Calculate Sparseness index for each cell
        for v = 1:numcell
           for j = 1:nFreq
               mFR_ON(n,j) = mean(FR_LASER_ALL(j,v,241:290));
               mFR_OFF(n,j) = mean(FR_NOLASER_ALL(j,v,241:290));
           end

            LASERSidxON(n, v) = Sparseness(mFR_ON(n,:),nFreq);  
            LASERSidxOFF(n, v) = Sparseness(mFR_OFF(n,:),nFreq);  
        end

    end

    descr{1} = {'Laser Onset:'; '-100 ms'};
    descr{2} = {'Laser Onset:'; '-20 ms'};
    descr{3} = {'Laser Onset:'; '+8 ms'};
    figure; 
    for n = 1:3
        subplot(3,2,(2*n)-1);
        scatter(LASERSidxOFF(n,:), LASERSidxON(n,:), 25,'filled')
        hold on; line([0 1], [0 1],'Color','k','LineStyle','--')
        hold off; axis square
        xlabel('Sparseness OFF')
        ylabel('Sparseness ON')
        set(gca,'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1])
        [p_LASERsparse(n), h_LASERsparse(n)] = signrank(LASERSidxOFF(n,:), LASERSidxON(n,:));
        text(-0.6,0.5, descr{n});
        title(['p = ' num2str(p_LASERsparse(n))]);
        
        subplot(3,2,2*(n));
        bar([nanmean(LASERSidxOFF(n,:)) nanmean(LASERSidxON(n,:))],0.5,'EdgeColor','none');
        hold on; errorbar([nanmean(LASERSidxOFF(n,:)) nanmean(LASERSidxON(n,:))],[nanstd(LASERSidxOFF(n,:))...
            ./sqrt(length(LASERSidxOFF(n,:))) nanstd(LASERSidxON(n,:))./sqrt(length(LASERSidxON(n,:)))],'k','LineStyle','none','LineWidth',2)
        box off
        ylabel('Sparseness')
        title(['N = ' num2str(numcell)]);
        set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out','YLim',[0 1]); axis square
    end


    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperUnits','points');
    set(gcf,'PaperSize',[800 1600]);
    set(gcf,'Position',[0 0 800 1600]);
    suptitle(['Feedback Affected Cells: ' TITLES{ii}])
    %cd(FigOut)
    %print(['sparsenessAFF_' num2str(ii)],'-dpdf','-r400')

end

%%
%%%%%%%%%%%%% III. LINEAR FITS %%%%%%%%%%%%%

onsets = {'-100 ms','-20 ms', '+8 ms'};

fitted = LinearFits.magDOWN;
for n = 1:3
    meanY = mean(fitted{n},1);
    semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
    subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
    hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
    line([0 1],[0 1],'Color','k','linestyle','--');
    plot([0:0.05:1],meanY + semY,'--r');
    plot([0:0.05:1],meanY - semY,'--r');
    hold off;
    box off;
    set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
    title(['Laser Onset: ' onsets{n}]);
end

suptitle('Feedback Affected Cells: Response Mag DECREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOut)
print('LinFit_magDOWN','-dpdf','-r400')


fitted = LinearFits.magUP;
for n = 1:3
    meanY = mean(fitted{n},1);
    semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
    subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
    hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
    line([0 1],[0 1],'Color','k','linestyle','--');
    plot([0:0.05:1],meanY + semY,'--r');
    plot([0:0.05:1],meanY - semY,'--r');
    hold off;
    box off;
    set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
    title(['Laser Onset: ' onsets{n}]);
end

suptitle('Feedback Affected Cells: Response Mag INCREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOut)
print('LinFit_magUP','-dpdf','-r400')


fitted = LinearFits.spontDOWN;
for n = 1:3
    meanY = mean(fitted{n},1);
    semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
    subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
    hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
    line([0 1],[0 1],'Color','k','linestyle','--');
    plot([0:0.05:1],meanY + semY,'--r');
    plot([0:0.05:1],meanY - semY,'--r');
    hold off;
    box off;
    set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
    title(['Laser Onset: ' onsets{n}]);
end

suptitle('Feedback Affected Cells: Response Spont DECREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOut)
print('LinFit_spontDOWN','-dpdf','-r400')


fitted = LinearFits.spontUP;
for n = 1:3
    meanY = mean(fitted{n},1);
    semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
    subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
    hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
    line([0 1],[0 1],'Color','k','linestyle','--');
    plot([0:0.05:1],meanY + semY,'--r');
    plot([0:0.05:1],meanY - semY,'--r');
    hold off;
    box off;
    set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
    title(['Laser Onset: ' onsets{n}]);
end

suptitle('Feedback Affected Cells: Response Spont INCREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOut)
print('LinFit_spontUP','-dpdf','-r400')


%% ************************************************************************
%  *****                         4. STRF FIGURES                      *****
%  ************************************************************************
