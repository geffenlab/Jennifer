clear all
FigOut = 'C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\Manuscripts\Feedback\Figures';
FileOut = 'C:\Users\Jennifer\Documents\MATLAB\TuningCurveAnalysis\';
tunedata = 'Feedback ChR2 in IC.mat';  %Calculated from TuningCurveAnalysis2.m
cd(FileOut)
load(tunedata);
opsin = 'chr2';

for z = 1:length(GOODCELL)
    len(z) = length(GOODCELL{z}(1,:));
end
maxlen = max(len);
tempcell = cell(1,length(GOODCELL));
for z = 1:length(GOODCELL)
    
    tempcell{z}(:,1:len(z)) = GOODCELL{z};
    if len(z) ~= maxlen
        tempcell{z}(:,len(z)+1:maxlen) = ' ';
    end
end
GOODCELLall = vertcat(tempcell{:});

% Use all cells (multiunits + single) to look for laser affect
for z = 1:length(allCELL)
    len(z) = length(allCELL{z}(1,:));
end
maxlen = max(len);
tempcell = cell(1,length(allCELL));
for z = 1:length(allCELL)
    
    tempcell{z}(:,1:len(z)) = allCELL{z};
    if len(z) ~= maxlen
        tempcell{z}(:,len(z)+1:maxlen) = ' ';
    end
end
temp = vertcat(tempcell{:});
Qcell = temp(horzcat(CellQ{:}) > 0 & horzcat(CellQ{:}) < 5,:);

fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k','defaultAxesFontSize',18);
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})


%% ************************************************************************
%  *****                       1. SILENCE + LASER                     *****
%  ************************************************************************


%%%%%%%%%%%%% I. EXAMPLE TRACE %%%%%%%%%%%%%
% Plot Raster and smoothed FR for both 25 ms and 250 ms stimuli (these were
% the two that significantly affect FR)

EX = {'M3222-2-9-1-1'; 'M3222-3-8-2-1'}; %Example units for ChR2
%EX = {}; %Example units for ArchT
for ex = 1:length(EX)
    cd('D:\Spikes')
    load(['R452-' EX{ex} '.mat']); %25 ms stim
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
    
    load(['R453-' EX{ex} '.mat']); %250 ms stim
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
    suptitle(EX{ex});
    
    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperUnits','points');
    set(gcf,'PaperSize',[1600 800]);
    set(gcf,'Position',[0 0 1600 800]);
    set(gcf,'renderer','painters');
    pause(1)
    cd(FigOut)
    print(['SilenceExample' num2str(ex) '_' opsin],'-dpdf','-r400')

end

%%%%%%%%%%%%% II. PLOT LASER ON VS LASER OFF RESPONSE %%%%%%%%%%%%%
cd('D:\Spikes')
Stim = {'R450';'R451';'R452';'R453'};
onset = 100; %Laser onset (ms)
dur = [1 5 25 250]; %Laser duration (ms) 
SilenceON = NaN(size(Qcell,1),length(Stim));
SilenceOFF = NaN(size(Qcell,1),length(Stim));
SilenceVAR = NaN(size(Qcell,1),length(Stim));

for u = 1:size(Qcell,1)
    
    for i = 1:length(Stim)
        q = find(Qcell(u,:) == '.');
        if exist([Stim{i} '-' Qcell(u,7:q - 10) '.mat'])
            load([Stim{i} '-' Qcell(u,7:q - 10) '.mat']);
        end
        if ~isempty(SpikeData)
            smoothedFR = smoothFRx4(SpikeData(3,:),nRep,0.001,[0 1],5);
            FRmeanON = mean(smoothedFR(onset:onset+dur(i)+20));
            FRmeanOFF = mean(smoothedFR(1:onset-1));
            FRsdevOFF = std(smoothedFR(750:end));
            if FRmeanOFF > 0.5 % Low spontaneous firing rates can lead to outliers
                SilenceON(u,i) = FRmeanON;
                SilenceOFF(u,i) = FRmeanOFF;
                SilenceVAR(u,i) = FRsdevOFF;
            end            
        end
    end
end
SilenceDEV = (SilenceON - SilenceOFF)./SilenceVAR;

figure;
edges = [-5:0.5:5];
%colour = [0.2 0.2 0.2; 0.4 0.4 0.4; 0.6 0.6 0.6; 0.8 0.8 0.8];
StimTitle = {'1 ms';'5 ms';'25 ms'; '250 ms'};
for i = 1:length(Stim)
    subplot(2,2,i);
    c = histc(SilenceDEV(:,i),edges);
    silence_mean(i) = nanmean(SilenceDEV(:,i));
    silence_stdev = nanstd(SilenceDEV(:,i));
    [h_silence(i),p_silence(i)] = kstest((SilenceDEV(:,i)-silence_mean(i))/silence_stdev);
    silence_median(i) = nanmedian(SilenceDEV(:,i));
    bar(edges,c,'EdgeColor','none')
    line([silence_mean(i) silence_mean(i)],[0 max(c)],'color','k','linestyle','--')
    line([silence_median(i) silence_median(i)],[0 max(c)],'color','r','linestyle','--')
    set(gca, 'XTick',edges(1:2:end),'XTickLabels',edges(1:2:end),'TickDir','out'); box off;
    xlabel('STDEVs above baseline')
    ylabel('cell count')
    title(StimTitle{i});
end
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 1000]);
set(gcf,'Position',[0 0 1600 1000]);
suptitle('Silence + Laser distributions')
%cd(FileOut)
%save(TITLE,'h_silence','p_silence','silence_mean','silence_median','-append')

cd(FigOut)
print(['SilenceDistrib_' opsin],'-dpdf','-r400')


%% ************************************************************************
%  *****                       2. CLICK FIGURES                       *****
%  ************************************************************************

ClickOFF = [0 2 4]; %Onset of clicks in laser off condition
ClickON = [1 3 5]; %Onset of clicks in laser on condition
ClickDUR = 0.05; %Window (in seconds) for single click duration
LaserON = [0.75 2.75 4.75]; %Laser onset (seconds)
LaserDUR = 1; %Laser duration (seconds)

%%%%%%%%%%%%% I. EXAMPLE TRACE %%%%%%%%%%%%%
cd('D:\Spikes')
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
set(gcf,'Renderer','painters');
cd(FigOut)
print(['ClickExample1_' opsin],'-dpdf','-r400')


%%%%%%%%%%%%% II. Laser ON vs Laser OFF RESPONSE %%%%%%%%%%%%%
cd('D:\Spikes')
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
[p_clickspont, h_clickspont] = signrank(spontOFF, spontON);
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
        ./sqrt(sum(isnan(ClickOFF))) nanstd(ClickON)./sqrt(sum(~isnan(ClickON)))],'k','LineStyle','none','LineWidth',2)
box off; ylabel('Firing Rate (Hz)'); title(['N = ' num2str(sum(~isnan(ClickOFF)))]);
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
suptitle('Click Response')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1000 600]);
set(gcf,'Position',[0 0 1000 600]);
%cd(FigOut)
%print(['ClickONOFF_' opsin],'-dpdf','-r400')

figure; subplot(1,2,1);
scatter(spontOFF,spontON,25,'filled','MarkerFaceColor',[0.2 0.2 0.2])
max1 = nanmax(spontOFF(:));
max2 = nanmax(spontON(:));
line([0 nanmax(max1,max2)],[0 nanmax(max1,max2)],'linestyle','--','color','k');
axis square; axis tight
set(gca,'TickDir','out');
xlabel('Firing Rate OFF (Hz)'); ylabel('Firing Rate ON (Hz)'); title(['p = ' num2str(p_clickspont)]);
subplot(1,2,2);
 bar([nanmean(spontOFF) nanmean(spontON)],0.5,'EdgeColor','none','FaceColor',[0.2 0.2 0.2]);
     hold on; errorbar([nanmean(spontOFF) nanmean(spontON)],[nanstd(spontOFF)...
        ./sqrt(sum(isnan(spontOFF))) nanstd(spontON)./sqrt(sum(~isnan(spontON)))],'k','LineStyle','none','LineWidth',2)
box off; ylabel('Firing Rate (Hz)'); title(['N = ' num2str(sum(~isnan(spontOFF)))]);
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
suptitle('Click Response (spontaneous)')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1000 600]);
set(gcf,'Position',[0 0 1000 600]);
cd(FigOut)
print(['spontClickONOFF_' opsin],'-dpdf','-r400')

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

%cd(FigOut)
%print(['ClickDELTA_' opsin],'-dpdf','-r400')


%% ************************************************************************
%  *****                       3. TUNING FIGURES                      *****
%  ************************************************************************

%%%%%%%%%%%%% I. EXAMPLE TCs AND TIMECOURSES %%%%%%%%%%%%%

tLaserON = [0.4 0.48 0.508];
tLaserOFF = [0.65 0.73 0.758];
onsets = {'-100 ms','-20 ms', '+8 ms'};
EX = {'M3222-1-5-3-1';'M3222-1-7-3-1'; 'M3223-1-8-3-1'; 'M3226-3-6-1-1'};
%EX = {'M3222-1-7-3-1'};
for ex = 2%1:length(EX)
    for i = 1:length(tLaserON)
        cd('D:\Spikes')
        load(['TC003-' EX{ex} '_laser-0' num2str(i)]);
        laserFR = vertcat(FR_LASER{:,i});
        nolaserFR = vertcat(FR_NOLASER{:,i});
        
        %Plot timecourses
        subplot(2, 3, i);
        %tt1 = find(TCon.t>=.2, 1, 'first');
        %tt5 = find(TCon.t>=.8, 1, 'first');
        %fr1 = SmoothGaus(mean(TCon.FRmat, 1), 3);
        %fr2 = SmoothGaus(mean(TCoff.FRmat, 1), 3);
        %plot(TCoff.t(tt1:tt5), fr1(tt1:tt5), 'k','linewidth',2);
        %hold on;
        %plot(TCoff.t(tt1:tt5), fr2(tt1:tt5), 'r','linewidth',2);
        plot([0.26:0.001:0.74-0.001],nolaserFR(23,:),'k','linewidth',2)
        hold on;
        plot([0.26:0.001:0.74-0.001],laserFR(23,:),'r','linewidth',2)
        ylabel('Firing rate (Hz)'); xlabel('Time (s)');
        line([tLaserON(i) tLaserOFF(i); tLaserON(i) tLaserOFF(i)], [0 0; max(laserFR(23,:)) max(laserFR(23,:))], 'color', 'm','linestyle','--');
        line([0.5 0.55; 0.5 0.55], [0 0; max(laserFR(23,:)) max(laserFR(23,:))], 'color', 'k','linestyle','--');
        box off; axis tight;
        hold off;
        set(gca,'TickDir','out','XTick',[0.3 0.5 0.7],'XLim',[0.24 0.76])
        title(['Laser Onset: ' onsets(i)]);

        Ton = SmoothGaus(TCon.TCmat{1}, 3);
        Toff = SmoothGaus(TCoff.TCmat{1}, 3);

        %Plot tuning curves
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
print(['TuneExample' num2str(ex) '_' opsin],'-dpdf','-r400')
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

for ii = 2%1:4
    h1 = cellidx{ii};
    numcell = size(h1,1);
    LASERSidxON = nan(3,numcell);
    LASERSidxOFF = nan(3,numcell);
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
    %suptitle(['Feedback Affected Cells: ' TITLES{ii}])
    %cd(FigOut)
    %print(['SparsenessAffected' num2str(ii) '_' opsin],'-dpdf','-r400')

end

%%
%%%%%%%%%%%%% III. LINEAR FITS %%%%%%%%%%%%%

onsets = {'-100 ms','-20 ms', '+8 ms'};

% % fitted = LinearFits.magDOWN;
% % for n = 1:3
% %     meanY = mean(fitted{n},1);
% %     semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
% %     subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
% %     hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
% %     line([0 1],[0 1],'Color','k','linestyle','--');
% %     plot([0:0.05:1],meanY + semY,'--r');
% %     plot([0:0.05:1],meanY - semY,'--r');
% %     hold off;
% %     box off;
% %     set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
% %     title(['Laser Onset: ' onsets{n}]);
% % end
% % 
% % suptitle('Feedback Affected Cells: Response Mag DECREASE')
% % set(gcf,'PaperPositionMode','auto');         
% % set(gcf,'PaperOrientation','landscape');
% % set(gcf,'PaperUnits','points');
% % set(gcf,'PaperSize',[1600 500]);
% % set(gcf,'Position',[0 0 1600 500]);
% % cd(FigOut)
% % print(['LinFit_magDOWN_' opsin],'-dpdf','-r400')
% % 
% % 
% % fitted = LinearFits.magUP;
% % for n = 1:3
% %     meanY = mean(fitted{n},1);
% %     semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
% %     subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
% %     hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
% %     line([0 1],[0 1],'Color','k','linestyle','--');
% %     plot([0:0.05:1],meanY + semY,'--r');
% %     plot([0:0.05:1],meanY - semY,'--r');
% %     hold off;
% %     box off;
% %     set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
% %     title(['Laser Onset: ' onsets{n}]);
% % end
% % 
% % suptitle('Feedback Affected Cells: Response Mag INCREASE')
% % set(gcf,'PaperPositionMode','auto');         
% % set(gcf,'PaperOrientation','landscape');
% % set(gcf,'PaperUnits','points');
% % set(gcf,'PaperSize',[1600 500]);
% % set(gcf,'Position',[0 0 1600 500]);
% % cd(FigOut)
% % print(['LinFit_magUP_' opsin],'-dpdf','-r400')
% % 
% % 
% % fitted = LinearFits.spontDOWN;
% % for n = 1:3
% %     meanY = mean(fitted{n},1);
% %     semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
% %     subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
% %     hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
% %     line([0 1],[0 1],'Color','k','linestyle','--');
% %     plot([0:0.05:1],meanY + semY,'--r');
% %     plot([0:0.05:1],meanY - semY,'--r');
% %     hold off;
% %     box off;
% %     set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
% %     title(['Laser Onset: ' onsets{n}]);
% % end
% % 
% % suptitle('Feedback Affected Cells: Response Spont DECREASE')
% % set(gcf,'PaperPositionMode','auto');         
% % set(gcf,'PaperOrientation','landscape');
% % set(gcf,'PaperUnits','points');
% % set(gcf,'PaperSize',[1600 500]);
% % set(gcf,'Position',[0 0 1600 500]);
% % cd(FigOut)
% % print(['LinFit_spontDOWN_' opsin],'-dpdf','-r400')
% % 
% % 
% % fitted = LinearFits.spontUP;
% % for n = 1:3
% %     meanY = mean(fitted{n},1);
% %     semY = std(fitted{n},[],1)./sqrt(size(fitted{n},1));
% %     subplot(1,3,n); plot([0:0.05:1],fitted{n}','Color', [0.8 0.8 0.8],'linewidth',2)
% %     hold on; plot([0:0.05:1],meanY,'Color','r','linewidth',2); 
% %     line([0 1],[0 1],'Color','k','linestyle','--');
% %     plot([0:0.05:1],meanY + semY,'--r');
% %     plot([0:0.05:1],meanY - semY,'--r');
% %     hold off;
% %     box off;
% %     set(gca,'TickDir','out'); xlabel('Normalized FR OFF'); ylabel('Normalized FR ON');
% %     title(['Laser Onset: ' onsets{n}]);
% % end
% % 
% % suptitle('Feedback Affected Cells: Response Spont INCREASE')
% % set(gcf,'PaperPositionMode','auto');         
% % set(gcf,'PaperOrientation','landscape');
% % set(gcf,'PaperUnits','points');
% % set(gcf,'PaperSize',[1600 500]);
% % set(gcf,'Position',[0 0 1600 500]);
% % cd(FigOut)
% % print(['LinFit_spontUP_' opsin],'-dpdf','-r400')

% Scatter linear fit slopes
figure;
for i = 1:length(linearparams)
    subplot(2,2,i)
    scatter(ones(1,length(linearparams{i})),linearparams{i}(:,1),'filled','k')
    line([1.1 1.3],[prctile(linearparams{i}(:,1),50) prctile(linearparams{i}(:,1),50)],'Color','r','LineWidth',1)
    line([1.2 1.2], [prctile(linearparams{i}(:,1),25) prctile(linearparams{i}(:,1),75)],'Color','r','LineWidth',2)

    hold on; scatter(2*ones(1,length(linearparams{i})),linearparams{i}(:,2),'filled','k')
    line([2.1 2.3],[prctile(linearparams{i}(:,2),50) prctile(linearparams{i}(:,2),50)],'Color','r','LineWidth',1)
    line([2.2 2.2], [prctile(linearparams{i}(:,2),25) prctile(linearparams{i}(:,2),75)],'Color','r','LineWidth',2)

    hold on; scatter(3*ones(1,length(linearparams{i})),linearparams{i}(:,3),'filled','k')
    line([3.1 3.3],[prctile(linearparams{i}(:,3),50) prctile(linearparams{i}(:,3),50)],'Color','r','LineWidth',1)
    line([3.2 3.2], [prctile(linearparams{i}(:,3),25) prctile(linearparams{i}(:,3),75)],'Color','r','LineWidth',2)

    set(gca,'TickDir','out','XTick', [1 2 3],'XTickLabel',[-100, -20, 8],'XLim',[0.5 3.5],'YLim',[0 3])
    line([0.5 3.5],[1 1],'LineStyle','--','Color','k')
    ylabel('Slope coefficients'); xlabel('Laser onset (ms)')

end

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 800]);
set(gcf,'Position',[0 0 800 800]);
%cd(FigOut)
%print(['Slope_coeff_' opsin],'-dpdf','-r400')

figure;
for i = 1:length(linearparams)
    subplot(2,2,i)
    scatter(ones(1,length(linearparams{i})),linearparams{i}(:,4),'filled','k')
    line([1.1 1.3],[prctile(linearparams{i}(:,4),50) prctile(linearparams{i}(:,4),50)],'Color','r','LineWidth',1)
    line([1.2 1.2], [prctile(linearparams{i}(:,4),25) prctile(linearparams{i}(:,4),75)],'Color','r','LineWidth',2)

    hold on; scatter(2*ones(1,length(linearparams{i})),linearparams{i}(:,5),'filled','k')
    line([2.1 2.3],[prctile(linearparams{i}(:,5),50) prctile(linearparams{i}(:,5),50)],'Color','r','LineWidth',1)
    line([2.2 2.2], [prctile(linearparams{i}(:,5),25) prctile(linearparams{i}(:,5),75)],'Color','r','LineWidth',2)

    hold on; scatter(3*ones(1,length(linearparams{i})),linearparams{i}(:,6),'filled','k')
    line([3.1 3.3],[prctile(linearparams{i}(:,6),50) prctile(linearparams{i}(:,6),50)],'Color','r','LineWidth',1)
    line([3.2 3.2], [prctile(linearparams{i}(:,6),25) prctile(linearparams{i}(:,6),75)],'Color','r','LineWidth',2)

    set(gca,'TickDir','out','XTick', [1 2 3],'XTickLabel',[-100, -20, 8],'XLim',[0.5 3.5])
    line([0.5 3.5],[0 0],'LineStyle','--','Color','k')
    ylabel('y-intercept'); xlabel('Laser onset (ms)')

end

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 800]);
set(gcf,'Position',[0 0 800 800]);
%cd(FigOut)
%print(['yint_coeff_' opsin],'-dpdf','-r400')

%%
%%%%%%%%%%%%% IV. PREFERRED FREQUENCY %%%%%%%%%%%%%

onsets = {'-100 ms','-20 ms', '+8 ms'};
cellidx{1} = GOODCELLall(SigCellmagDOWN,:);
cellidx{2} = GOODCELLall(SigCellmagUP,:);
cellidx{3} = GOODCELLall(SigCellspontDOWN,:);
cellidx{4} = GOODCELLall(SigCellspontUP,:);

nfreq = 1;
amps = 1;
for stimidx = 1:length(onsets)
    for u = 1:length(cellidx)
        h1 = cellidx{u};
        for v = 1:size(h1,1)
            q = find(h1(v,:) == '.');
            load([h1(v,1:q - 2) num2str(stimidx) '.mat']);

            bestfreq.ON{u}(v) = calcPreferredFreqsTC(TCon, nfreq, amps);
            bestfreq.OFF{u}(v) = calcPreferredFreqsTC(TCoff, nfreq, amps);
        end
        f1 = figure(1);

        subplot(3,4,4*(stimidx-1) + u);
        scatter(bestfreq.OFF{u}/1000,bestfreq.ON{u}/1000,25,'filled','b')
        
        hold on; line([0 70], [0 70],'color','k', 'linestyle','--');
        axis square; axis tight;
        set(gca,'TickDir','out','XTick',[0:20:60],'YTick',[0:20:60])
        
        f2 = figure(2);

        subplot(3,4,4*(stimidx-1) + u);
        plotErrorBar([[bestfreq.OFF{u}/1000]' [bestfreq.ON{u}/1000]'],'sem',{'OFF'; 'ON'});
        set(gca,'TickDir','out')

        
        [p_BF(u,stimidx),h] = signrank(bestfreq.OFF{u}, bestfreq.ON{u});
    end

    
end
set(f1,'PaperPositionMode','auto');         
set(f1,'PaperOrientation','landscape');
set(f1,'PaperUnits','points');
set(f1,'PaperSize',[1200 800]);
set(f1,'Position',[0 0 1200 800]);

set(f2,'PaperPositionMode','auto');         
set(f2,'PaperOrientation','landscape');
set(f2,'PaperUnits','points');
set(f2,'PaperSize',[1200 800]);
set(f2,'Position',[0 0 1200 800]);
%cd(FigOut)
%print(f1, ['bestfreqSCAT_' opsin],'-dpdf','-r400')
%print(f2, ['bestfreqBAR_' opsin],'-dpdf','-r400')

%%
%%%%%%%%%%%%% V. SPONTANEOUS AND TONE-EVOKED RESPONSE MAGNITUDE %%%%%%%%%%%%%
onsets = [-100 -20 8];
for stimidx = 1:length(onsets)

spontIDX = 219:239; %Indices for spontaneous response
toneIDX = 240:315; %Indices for tone-evoked response
spontON = cell(1,length(MNum));
spontOFF = cell(1,length(MNum));
aveON = cell(1,length(MNum));
aveOFF = cell(1,length(MNum));
magON = cell(1,length(MNum));
magOFF = cell(1,length(MNum));


for u = 1:length(MNum)
    
    %Calculate spontaneous firing rate and tone-evoked response magnitude
    numcellB = length(IDX{u});
    spontON{u} = NaN(1,numcellB); spontOFF{u} = NaN(1,numcellB); magON{u} = NaN(1,numcellB);
    magOFF{u} = NaN(1,numcellB); aveON{u} = NaN(1,numcellB); aveOFF{u} = NaN(1,numcellB);
    for j = 1:numcellB   
        spontON{u}(j) = mean(FR_LASER{u,stimidx}(IDX{u}(j),spontIDX)); spontOFF{u}(j) = mean(FR_NOLASER{u,stimidx}(IDX{u}(j),spontIDX));
        aveON{u}(j) = mean(FR_LASER{u,stimidx}(IDX{u}(j),toneIDX)); aveOFF{u}(j) = mean(FR_NOLASER{u,stimidx}(IDX{u}(j),toneIDX));
        magON{u}(j) = aveON{u}(j) - spontON{u}(j); magOFF{u}(j) = aveOFF{u}(j) - spontOFF{u}(j);
    end
    
end

spontOFFall = [spontOFF{:}];
spontONall = [spontON{:}];
magOFFall = [magOFF{:}];
magONall = [magON{:}];
aveONall = [aveON{:}];
aveOFFall = [aveOFF{:}];

% % figure;
% % subplot(1,4,1); scatter(spontOFFall, spontONall,25,'filled')
% % LIM = max([spontOFFall, spontONall]);
% % %LIM = 50;
% % hold on; line([0 LIM], [0 LIM],'Color','k','LineStyle','--');
% % hold off;
% % set(gca,'TickDir','out'); axis square
% % set(gca,'xlim',[0 LIM],'ylim',[0 LIM])
% % xlabel('LASER OFF')
% % ylabel('LASER ON')
% % title('Spontaneous activity (Hz)')
% % 
% % subplot(1,4,2); bar([mean(spontOFFall) mean(spontONall)],0.5,'EdgeColor','none');
% % hold on; errorbar([mean(spontOFFall) mean(spontONall)],[std(spontOFFall)./sqrt(length(spontOFFall)) std(spontONall)./sqrt(length(spontONall))],'k','LineStyle','none','LineWidth',2)
% % box off
% % ylabel('Firing rate (Hz)')
% % set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
% % %set(gca,'YLim',[0 15])
% % [p_spont(stimidx), h_spont(stimidx)] = signrank(spontOFFall, spontONall);
% % 
% % subplot(1,4,3); scatter(magOFFall, magONall,25,'filled')
% % LIM = max([magOFFall, magONall]);
% % %LIM = 50;
% % hold on; line([0 LIM], [0 LIM], 'Color','k','LineStyle','--');
% % hold off;
% % set(gca,'TickDir','out'); axis square
% % set(gca,'xlim',[0 LIM],'ylim',[0 LIM])
% % xlabel('LASER OFF')
% % ylabel('LASER ON')
% % title('Tone-evoked response magnitude (Hz)')
% % 
% % subplot(1,4,4); bar([mean(magOFFall) mean(magONall)],0.5,'EdgeColor','none');
% % hold on; errorbar([mean(magOFFall) mean(magONall)],[std(magOFFall)./sqrt(length(magOFFall)) std(magONall)./sqrt(length(magONall))],'k','LineStyle','none','LineWidth',2)
% % box off
% % ylabel('Firing rate (Hz)')
% % set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
% % %set(gca,'YLim',[0 30])
% % [p_mag(stimidx), h_mag(stimidx)] = signrank(magOFFall, magONall);
% % 
% % suptitle({TITLE, ['Laser Onset: ' num2str(onsets(stimidx))]});

[p_ave(stimidx), h_ave(stimidx)] = signrank(aveOFFall, aveONall);
figure; subplot(1,2,1); scatter(aveOFFall,aveONall,25,'filled');
LIM = max([aveOFFall, aveONall]);
hold on; line([0 LIM], [0 LIM], 'Color','k','LineStyle','--');
xlabel('LASER OFF')
ylabel('LASER ON')
title(['p = ' num2str(p_ave(stimidx))])
hold off;
set(gca,'TickDir','out'); axis square
set(gca,'xlim',[0 LIM],'ylim',[0 LIM])
subplot(1,2,2); plotErrorBar([aveOFFall' aveONall'],'sem',{'OFF';'ON'})
box off
ylabel('Firing rate (Hz)')
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
title('Tone-evoked response (Hz)')
suptitle({TITLE, ['Laser Onset: ' num2str(onsets(stimidx))]});
%FileOut = 'C:\Users\Jennifer\Documents\MATLAB\TuningCurveAnalysis';
%cd(FileOut);
%save(TITLE,'h_mag','p_mag','h_spont','p_spont','p_ave','-append')
aveON = mean(aveOFFall)
aveOFF = mean(aveONall)
aveONsem = nanstd(aveONall)./sqrt(length(find(~isnan(aveONall))))
aveOFFsem = nanstd(aveOFFall)./sqrt(length(find(~isnan(aveOFFall))))

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 600]);
set(gcf,'Position',[0 0 1600 600]);

%cd(FigOut)
%print([TITLE '_' num2str(stimidx)],'-dpdf','-r400')
end                                                                                                                                                                         


%% ************************************************************************
%  *****                         4. STRF FIGURES                      *****
%  ************************************************************************

%%%%%%%%%%%%% I. EXAMPLE STRFs %%%%%%%%%%%%%
%EX = {'M3222-1-6-1-1';'M3222-1-7-3-1';'M3222-3-3-2-1';'M3222-3-6-1-1'}; %Example increasing units
EX = {'M3222-1-7-1-1';'M3222-4-2-2-1'}; %Example decreasing units
load('DRC001-01','params');

for ex = 1:length(EX);
    cd('D:\Spikes')
    load(['DRC001-' EX{ex} '.mat']);
    
    MM = max(max([STAon;STAoff3]));
    mm = min(min([STAon;STAoff3]));
    subplot(1,2,1); h = plotSTA([-0.1:0.005:0],params.freqs/1000,STAoff3,1,[mm,MM]);
    title('Laser OFF');
    
    
    subplot(1,2,2); h = plotSTA([-0.1:0.005:0],params.freqs/1000,STAon,1,[mm,MM]); 
    title('Laser ON');
    suptitle(EX{ex})
    colorbar
    
    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperUnits','points');
    set(gcf,'PaperSize',[1200 600]);
    set(gcf,'Position',[0 0 1200 600]);
    set(gcf,'Renderer','painters');
    
    cd(FigOut)
    print(['STRFexample' num2str(ex) '_' opsin],'-dpdf','-r400')    
end



%%
%%%%%%%%%%%%% II. PLOT STRF PARAMETERS LASER ON VS LASER OFF %%%%%%%%%%%%%
%Sorry for this ugly code

paramlabels = {'time width (s)','freq width (kHz)','peak time (s)','peak freq (kHz)','size (pixels)', 'mean value'};
%For cells that show increase in FR with laser ON
POSparamsON = [];
POSparamsOFF = [];
NEGparamsON = [];
NEGparamsOFF = [];
noClust = [0 0 0 NaN NaN 0 0]; %This will be the data if a cluster in laser OFF is not present in laser ON

cd('D:\Spikes')
for u = 1:size(strfcellUP,1)
    q = find(strfcellUP(u,:) == '.');
    load(['DRC001-' strfcellUP(u,7:q-10) '.mat']);
    %First need to match clusters in OFF and ON conditions
    clustOFF = STRFclustOFF.POS.params.data(:,1);
    if clustOFF ~=0
        tempOFFparams = STRFclustOFF.POS.params.data;
        tempONparams = repmat(noClust,length(clustOFF),1);
        clustON = STRFclustON.POS.params.data(:,1);
        counter = 0;
        for w = clustOFF'
            counter = counter + 1;
            a = find(clustON == w);
            if ~isempty(a)
                %tempONparams(counter,:) = [STRFclustON.POS.params.data(a,1:end-1) STRFclustON.POS.params.OFFmaskmean(find(STRFclustON.POS.params.OFFmaskmean(:,1) == w,1,'first'),2)] ;
                tempONparams(counter,:) = STRFclustON.POS.params.data(a,1:end);
            end

        end

        POSparamsOFF = [POSparamsOFF; tempOFFparams];
        POSparamsON = [POSparamsON; tempONparams];
    end
    
    clustOFF = STRFclustOFF.NEG.params.data(:,1);
    if clustOFF ~=0
        tempOFFparams = STRFclustOFF.NEG.params.data;
        tempONparams = repmat(noClust,length(clustOFF),1);
        clustON = STRFclustON.NEG.params.data(:,1);
        counter = 0;
        for w = clustOFF'
            counter = counter + 1;
            a = find(clustON == w);
            if ~isempty(a)
                %tempONparams(counter,:) = [STRFclustON.NEG.params.data(a,1:end-1) STRFclustON.NEG.params.OFFmaskmean(find(STRFclustON.NEG.params.OFFmaskmean(:,1) == w,1,'first'),2)] ;
                tempONparams(counter,:) = STRFclustON.NEG.params.data(a,1:end);
            end

        end

        NEGparamsOFF = [NEGparamsOFF; tempOFFparams];
        NEGparamsON = [NEGparamsON; tempONparams];
    end
end

idxPOSgone = find(isnan(POSparamsON(:,4)));
percentPOSgone = (length(idxPOSgone)./length(POSparamsON(:,4)))*100;
idxNEGgone = find(isnan(NEGparamsON(:,4)));
percentNEGgone = (length(idxNEGgone)./length(NEGparamsON(:,4)))*100;

POSparamsOFF2 = POSparamsOFF;
POSparamsOFF2(idxPOSgone,:) = [];
POSparamsON2 = POSparamsON;
POSparamsON2(idxPOSgone,:) = [];

NEGparamsOFF2 = NEGparamsOFF;
NEGparamsOFF2(idxNEGgone,:) = [];
NEGparamsON2 = NEGparamsON;
NEGparamsON2(idxNEGgone,:) = [];

figure;
for i = 1:6
    subplot(2,3,i); scatter(POSparamsOFF2(:,i+1),POSparamsON2(:,i+1),30,'filled','d','markerfacecolor','k')
    hold on; scatter(NEGparamsOFF2(:,i+1),NEGparamsON2(:,i+1),25,'filled','markerfacecolor',[0.4 0.4 0.4]);
    MM = nanmax([POSparamsOFF2(:,i+1);POSparamsON2(:,i+1);NEGparamsOFF2(:,i+1);NEGparamsON2(:,i+1)]);
    mm = nanmin([POSparamsOFF2(:,i+1);POSparamsON2(:,i+1);NEGparamsOFF2(:,i+1);NEGparamsON2(:,i+1)]);
    line([mm MM],[mm MM],'Color','k','linestyle','--');
    [p_POS(i),~] = signrank(POSparamsOFF2(:,i+1),POSparamsON2(:,i+1));
    [p_NEG(i),~] = signrank(NEGparamsOFF2(:,i+1),NEGparamsON2(:,i+1));
    axis tight; axis square; box off;
    set(gca,'TickDir','out');
    xlabel('Laser OFF'); ylabel('Laser ON');
    title(paramlabels{i});
    legend(['p POS = ' num2str(p_POS(i))],['p NEG = ' num2str(p_NEG(i))],'location','best')

end
suptitle('STRF parameters: units with INCREASE in FR')

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 800]);
set(gcf,'Position',[0 0 1600 800]);
%cd(FigOut)
%print(['STRFparamsUP_' opsin],'-dpdf','-r400')  

figure;
for i = 1:6
    subplot(2,6,i); plotErrorBar([POSparamsOFF2(:,i+1) POSparamsON2(:,i+1)],'sem',{'OFF'; 'ON'});
    title(paramlabels{i});
    subplot(2,6,i+6); plotErrorBar([NEGparamsOFF2(:,i+1) NEGparamsON2(:,i+1)],'sem',{'OFF'; 'ON'});
end
suptitle('STRF parameters: units with INCREASE in FR')

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[2000 800]);
set(gcf,'Position',[0 0 2000 800]);
%cd(FigOut)
%print(['STRFbarUP_' opsin],'-dpdf','-r400')  

figure; h1 = subplot(2,1,1);
pie([percentPOSgone, 100-percentPOSgone],{['gone (' num2str(percentPOSgone) ,'%)'],['stable (' num2str(100-percentPOSgone) ,'%)']});
colormap(h1,[1 1 1;0 0 0])
title('Positive lobes')
h2 = subplot(2,1,2);
pie([percentNEGgone, 100-percentNEGgone],{['gone (' num2str(percentNEGgone) ,'%)'],['stable (' num2str(100-percentNEGgone) ,'%)']});
colormap(h2,[1 1 1;0.4 0.4 0.4])
title('Negative lobes')
suptitle('INCREASE FR units')

set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[600 1600]);
set(gcf,'Position',[0 0 600 1600]);
%cd(FigOut)
%print(['STRFpieUP_' opsin],'-dpdf','-r400') 

POSparamsOFF = [];
POSparamsON = [];
NEGparamsOFF = [];
NEGparamsON = [];
noClust = [0 0 0 NaN NaN 0 0]; %This will be the data if a cluster in laser OFF is not present in laser ON

cd('D:\Spikes')
for u = 1:size(strfcellDOWN,1)
    q = find(strfcellDOWN(u,:) == '.');
    load(['DRC001-' strfcellDOWN(u,7:q-10) '.mat']);
    
    %First need to match clusters in OFF and ON conditions
    clustOFF = STRFclustOFF.POS.params.data(:,1);
    if clustOFF ~= 0
        tempOFFparams = STRFclustOFF.POS.params.data;
        tempONparams = repmat(noClust,length(clustOFF),1);
        clustON = STRFclustON.POS.params.data(:,1);
        counter = 0;
        for w = clustOFF'
            counter = counter + 1;
            a = find(clustON == w);
            if ~isempty(a)
                %tempONparams(counter,:) = [STRFclustON.POS.params.data(a,1:end-1) STRFclustON.POS.params.OFFmaskmean(find(STRFclustON.POS.params.OFFmaskmean(:,1) == w,1,'first'),2)] ;
                tempONparams(counter,:) = [STRFclustON.POS.params.data(a,1:end)];
            end

        end

        POSparamsOFF = [POSparamsOFF; tempOFFparams];
        POSparamsON = [POSparamsON; tempONparams];
    end


    clustOFF = STRFclustOFF.NEG.params.data(:,1);
    if clustOFF ~= 0
        tempOFFparams = STRFclustOFF.NEG.params.data;
        tempONparams = repmat(noClust,length(clustOFF),1);
        clustON = STRFclustON.NEG.params.data(:,1);
        counter = 0;
        for w = clustOFF'
            counter = counter + 1;
            a = find(clustON == w);
            if ~isempty(a)
                %tempONparams(counter,:) = [STRFclustON.NEG.params.data(a,1:end-1) STRFclustON.NEG.params.OFFmaskmean(find(STRFclustON.NEG.params.OFFmaskmean(:,1) == w,1,'first'),2)] ;
                tempONparams(counter,:) = [STRFclustON.NEG.params.data(a,1:end)];
            end

        end

        NEGparamsOFF = [NEGparamsOFF; tempOFFparams];
        NEGparamsON = [NEGparamsON; tempONparams];
    end
end


idxPOSgone = find(isnan(POSparamsON(:,4)));
percentPOSgone = (length(idxPOSgone)./length(POSparamsON(:,4)))*100;
idxNEGgone = find(isnan(NEGparamsON(:,4)));
percentNEGgone = (length(idxNEGgone)./length(NEGparamsON(:,4)))*100;

POSparamsOFF2 = POSparamsOFF;
POSparamsOFF2(idxPOSgone,:) = [];
POSparamsON2 = POSparamsON;
POSparamsON2(idxPOSgone,:) = [];

NEGparamsOFF2 = NEGparamsOFF;
NEGparamsOFF2(idxNEGgone,:) = [];
NEGparamsON2 = NEGparamsON;
NEGparamsON2(idxNEGgone,:) = [];


figure;
for i = 1:6
    subplot(2,3,i); scatter(POSparamsOFF2(:,i+1),POSparamsON2(:,i+1),30,'d','filled','markerfacecolor','k')
    hold on; scatter(NEGparamsOFF2(:,i+1),NEGparamsON2(:,i+1),25,'filled','markerfacecolor',[0.4 0.4 0.4]);
    MM = nanmax([POSparamsOFF2(:,i+1);POSparamsON2(:,i+1);NEGparamsOFF2(:,i+1);NEGparamsON2(:,i+1)]);
    mm = nanmin([POSparamsOFF2(:,i+1);POSparamsON2(:,i+1);NEGparamsOFF2(:,i+1);NEGparamsON2(:,i+1)]);
    [p_POS(i),~] = signrank(POSparamsOFF2(:,i+1),POSparamsON2(:,i+1));
    [p_NEG(i),~] = signrank(NEGparamsOFF2(:,i+1),NEGparamsON2(:,i+1));
    line([mm MM],[mm MM],'Color','k','linestyle','--');
    axis tight; axis square; box off;
    set(gca,'TickDir','out');
    xlabel('Laser OFF'); ylabel('Laser ON');
    %title({paramlabels{i};['p = ' num2str(p_POS(i))]});
    title(paramlabels{i})
    legend(['p POS = ' num2str(p_POS(i))],['p NEG = ' num2str(p_NEG(i))],'location','best')
end
suptitle('STRF parameters: units with DECREASE in FR')


set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 1000]);
set(gcf,'Position',[0 0 1600 1000]);
%cd(FigOut)
%print(['STRFparamsDOWN_' opsin],'-dpdf','-r400') 

figure;
for i = 1:6
    subplot(2,6,i); plotErrorBar([POSparamsOFF2(:,i+1) POSparamsON2(:,i+1)],'sem',{'OFF'; 'ON'});
    title(paramlabels{i});
    subplot(2,6,i+6); plotErrorBar([NEGparamsOFF2(:,i+1) NEGparamsON2(:,i+1)],'sem',{'OFF'; 'ON'});
end
suptitle('STRF parameters: units with DECREASE in FR')

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[2000 800]);
set(gcf,'Position',[0 0 2000 800]);
%cd(FigOut)
%print(['STRFbarDOWN_' opsin],'-dpdf','-r400')  

figure; h1 = subplot(2,1,1);
pie([percentPOSgone, 100-percentPOSgone],{['gone (' num2str(percentPOSgone) ,'%)'],['stable (' num2str(100-percentPOSgone) ,'%)']});
colormap(h1,[1 1 1;0 0 0])
title('Positive lobes')
h2 = subplot(2,1,2);
pie([percentNEGgone, 100-percentNEGgone],{['gone (' num2str(percentNEGgone) ,'%)'],['stable (' num2str(100-percentNEGgone) ,'%)']});
colormap(h2,[1 1 1;0.4 0.4 0.4])
title('Negative lobes')
suptitle('DECREASE FR units')

set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[600 1600]);
set(gcf,'Position',[0 0 600 1600]);
%cd(FigOut)
%print(['STRFpieDOWN_' opsin],'-dpdf','-r400') 

%% Compare sparseness to change in lobe size

paramlabels = {'time width (s)','freq width (kHz)','peak time (s)','peak freq (kHz)','size (pixels)', 'mean value'};
%For cells that show increase in FR with laser ON
POSparamsON = [];
POSparamsOFF = [];
noClust = [NaN 0]; %This will be the data if a cluster in laser OFF is not present in laser ON
paramID = 6; %Parameter 6 is size

Win = 0.24; %Window (in seconds) around start of tone to look for spikes
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
LASERSidxON = cell(1,3);
LASERSidxOFF = cell(1,3);

onsets = [-100 -20 8];
cd('D:\Spikes')
for u = 1:size(strfcellUP,1)
    q = find(strfcellUP(u,:) == '.');
    load(['DRC001-' strfcellUP(u,7:q-10) '.mat']);
    
    %First need to match clusters in OFF and ON conditions
    clustOFF = STRFclustOFF.POS.params.data(:,1);
    if clustOFF ~=0
        tempOFFparams = STRFclustOFF.POS.params.data(:,paramID);
        tempONparams = repmat(noClust,length(clustOFF),1);
        clustON = STRFclustON.POS.params.data(:,1);
        counter = 0;
        for w = clustOFF'
            counter = counter + 1;
            a = find(clustON == w);
            if ~isempty(a)
                tempONparams(counter,:) = STRFclustON.POS.params.data(a,[4 paramID]);
            end

        end
        POSparamsOFF = [POSparamsOFF; tempOFFparams];
        POSparamsON = [POSparamsON; tempONparams];

 
        for n = 1:3


            SpkTime_LaserAll = cell(1,nFreq);
            SpkTime_NoLaserAll = cell(1,nFreq);

            load(['D:\Spikes\M' strfcellUP(u,8:11) '\SpikeMat\TC003_LEFT-0' num2str(n) '-' strfcellUP(u,7:q-10) '.mat']); 

            SpkTime_NoLaser = SpikeTime(StimOrder_NoLaser,SpikeData,nRep, Win);
            SpkTime_Laser = SpikeTime(StimOrder_Laser,SpikeData,nRep, Win);

            for j = 1:nFreq
                SpkTime_LaserAll{j} = SpkTime_Laser(j,amps,:);
                SpkTime_NoLaserAll{j} = SpkTime_NoLaser(j,amps,:);
            end

            %Calculate and smooth firing rate 
            SpikeDataLaserAll = cell(1,nFreq);
            SpikeDataNoLaserAll = cell(1,nFreq);
            clear LASER_SPIKE_ALL NOLASER_SPIKE_ALL
            %For laser trials
          for j = 1:nFreq
                    for w = 1:numel(SpkTime_LaserAll{j})
                        LASER_SPIKE_ALL{j}{w} = [sort(SpkTime_LaserAll{j}{w}); w*ones(1,length(SpkTime_LaserAll{j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                        SpikeDataLaserAll{j} = [SpikeDataLaserAll{j} LASER_SPIKE_ALL{j}{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
                    end
                    FR_LASER_ALL(j,:) = smoothFRx4(SpikeDataLaserAll{j},numel(LASER_SPIKE_ALL{j})*stimInfo.repeats,0.001,[-Win Win],5);

                    %For no laser trials
                    for w = 1:numel(SpkTime_NoLaserAll{j})
                        NOLASER_SPIKE_ALL{j}{w} = [sort(SpkTime_NoLaserAll{j}{w}); w*ones(1,length(SpkTime_NoLaserAll{j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                        SpikeDataNoLaserAll{j} = [SpikeDataNoLaserAll{j} NOLASER_SPIKE_ALL{j}{w}];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
                    end
                    FR_NOLASER_ALL(j,:) = smoothFRx4(SpikeDataNoLaserAll{j},numel(NOLASER_SPIKE_ALL{j})*stimInfo.repeats,0.001,[-Win Win],5);
           end


            %Calculate Sparseness index for each cell

           for j = 1:nFreq
               mFR_ON(n,j) = mean(FR_LASER_ALL(j,241:290));
               mFR_OFF(n,j) = mean(FR_NOLASER_ALL(j,241:290));
           end

           for i = clustOFF'
                LASERSidxON{n} = [LASERSidxON{n}; Sparseness(mFR_ON(n,:),nFreq)];  
                LASERSidxOFF{n} = [LASERSidxOFF{n}; Sparseness(mFR_OFF(n,:),nFreq)];  
           end

        end
    end
end

idxPOSgone = find(isnan(POSparamsON(:,1)));

POSparamsOFF2 = POSparamsOFF;
POSparamsOFF2(idxPOSgone,:) = [];
POSparamsON2 = POSparamsON;
POSparamsON2(idxPOSgone,:) = [];
STRFdiff = POSparamsON2(:,2) - POSparamsOFF2;

for n = 1:3
    LASERSidxON{n}(idxPOSgone,:) = [];
    LASERSidxOFF{n}(idxPOSgone,:) = [];
    SPARSEdiff{n} = LASERSidxON{n} - LASERSidxOFF{n};
end


figure; 
for n = 1:3
    subplot(1,3,n);
    scatter(STRFdiff,SPARSEdiff{n},'filled')
    hold on; line([0 0],[-0.5 0.5],'color','k','linestyle','--');
    line([min(STRFdiff) max(STRFdiff)],[0 0],'color','k','linestyle','--');
    set(gca,'TickDir','out','YLim',[-0.5 0.5])
    xlabel('delta STRF size (ON-OFF)');
    ylabel('delta Sparseness (ON - OFF)');
    title(['Laser Onset: ' num2str(onsets(n)) ' ms']);
    axis square
    axis tight
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 600]);
set(gcf,'Position',[0 0 1600 600]);
%cd(FigOut)
%print(['STRFvSparseness_' opsin],'-dpdf','-r400') 
