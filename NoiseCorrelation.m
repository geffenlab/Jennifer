%Calculate noise correlations
FileOut = 'C:\Users\Jennifer\Documents\MATLAB\TuningCurveAnalysis\';
FigOut = 'C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\Manuscripts\Feedback\Figures';
tunedata = 'Feedback ChR2 in IC.mat';  %Calculated from TuningCurveAnalysis2.m
cd(FileOut)
load(tunedata);

%% ***********************************************************
%  *****            1. CALCULATE SPIKE COUNTS            *****
%  ***********************************************************

%Load stimulus parameters
load('D:\Code\TuningCurve\TC003_170815LEFT_stimInfo');
trialdur = (stimInfo.tDur+stimInfo.ITI)/1000;

WinSpont = [-0.02 0];
WinTone = [0 0.075];
SpkCount_NoLaser_Spont = cell(1,length(GOODCELL)); SpkCount_Laser_Spont = cell(1,length(GOODCELL)); 
SpkCount_NoLaser_Tone = cell(1,length(GOODCELL)); SpkCount_Laser_Tone = cell(1,length(GOODCELL));


for u = 1:length(GOODCELL)
    SessionNums = str2num(GOODCELL{u}(:,13)); %Extract session ID for each neuron
    recSessions = unique(SessionNums); %Identify number of distinct recording sessions for each mouse
    nSessions = length(recSessions);
    SpkCount_NoLaser_Spont{u} = cell(1,nSessions); SpkCount_Laser_Spont{u} = cell(1,nSessions); 
    SpkCount_NoLaser_Tone{u} = cell(1,nSessions); SpkCount_Laser_Tone{u} = cell(1,nSessions);
    for v = 1:nSessions
        tempcells = GOODCELL{u}(find(SessionNums == v),:);
        nCells = size(tempcells,1);
        
        for i = 1:nCells
            match = strfind(tempcells(i,:), '.');
            q = match(end);
            load(['D:\Spikes\' tempcells(i,7:11) '\SpikeMat\TC003_LEFT-01' tempcells(i,6:q-10) '.mat']); 
            
            dur = trialdur*length(stimInfo.order)*nRep;
            Time = 0:trialdur:dur; %Each repetition is 400 seconds long, each trial is 500ms long
            Time = Time(1:end-1);
            StimOrder_Laser = [Time(2:2:end)];
            StimOrder_NoLaser = [Time(1:2:end)];
            
            %Extract spike counts
            SpkCount_NoLaser_Spont{u}{v}(:,i) = getSpikeCount_Trial(StimOrder_NoLaser,SpikeData, WinSpont);
            SpkCount_Laser_Spont{u}{v}(:,i) = getSpikeCount_Trial(StimOrder_Laser,SpikeData, WinSpont);
            SpkCount_NoLaser_Tone{u}{v}(:,i) = getSpikeCount_Trial(StimOrder_NoLaser,SpikeData, WinTone);
            SpkCount_Laser_Tone{u}{v}(:,i) = getSpikeCount_Trial(StimOrder_Laser,SpikeData, WinTone);
        end
        
    end
end

%% ***********************************************************
%  *****         2. CALCULATE NOISE CORRELATIONS         *****
%  ***********************************************************
corr_NoLaser_Spont = [];
corr_Laser_Spont = [];
corr_NoLaser_Tone = [];
corr_Laser_Tone = [];
for u = 1:length(SpkCount_NoLaser_Spont)
    for v = 1:length(SpkCount_NoLaser_Spont{u})
        if u == 1
            idx = v;
            corr_NoLaser_Spont{v} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
            corr_Laser_Spont{v} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
            corr_NoLaser_Tone{v} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
            corr_Laser_Tone{v} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
        else 
            idx = 1 + length(corr_NoLaser_Spont);
            corr_NoLaser_Spont{1+length(corr_NoLaser_Spont)} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
            corr_Laser_Spont{1+length(corr_Laser_Spont)} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
            corr_NoLaser_Tone{1+length(corr_NoLaser_Tone)} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
            corr_Laser_Tone{1+length(corr_Laser_Tone)} = zeros(size(SpkCount_NoLaser_Spont{u}{v},2),size(SpkCount_NoLaser_Spont{u}{v},2));
        end
        
        for i = 1:size(SpkCount_NoLaser_Spont{u}{v},2)            
            for j = i:size(SpkCount_NoLaser_Spont{u}{v},2)
                corr = corrcoef(SpkCount_NoLaser_Spont{u}{v}(:,i),SpkCount_NoLaser_Spont{u}{v}(:,j));
                corr_NoLaser_Spont{idx}(i,j) = corr(2);
                
                corr = corrcoef(SpkCount_Laser_Spont{u}{v}(:,i),SpkCount_Laser_Spont{u}{v}(:,j));
                corr_Laser_Spont{idx}(i,j) = corr(2);
                
                corr = corrcoef(SpkCount_NoLaser_Tone{u}{v}(:,i),SpkCount_NoLaser_Tone{u}{v}(:,j));
                corr_NoLaser_Tone{idx}(i,j) = corr(2);
                
                corr = corrcoef(SpkCount_Laser_Tone{u}{v}(:,i),SpkCount_Laser_Tone{u}{v}(:,j));
                corr_Laser_Tone{idx}(i,j) = corr(2);
            end
        end
        %Turn upper triangular matrix into full matrix
        corr_NoLaser_Spont{idx} = corr_NoLaser_Spont{idx} + corr_NoLaser_Spont{idx}' - diag(diag(corr_NoLaser_Spont{idx}));
        corr_Laser_Spont{idx} = corr_Laser_Spont{idx} + corr_Laser_Spont{idx}' - diag(diag(corr_Laser_Spont{idx}));
        corr_NoLaser_Tone{idx} = corr_NoLaser_Tone{idx} + corr_NoLaser_Tone{idx}' - diag(diag(corr_NoLaser_Tone{idx}));
        corr_Laser_Tone{idx} = corr_Laser_Tone{idx} + corr_Laser_Tone{idx}' - diag(diag(corr_Laser_Tone{idx}));

    end
end

%% ***********************************************************
%  *****           3. PLOT NOISE CORRELATIONS            *****
%  ***********************************************************
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k','defaultAxesFontSize',18);
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

nanjet = [ 1,1,1; jet ];

for v = 1:length(corr_NoLaser_Spont)
    figure;
    subplot(2,3,1); imagesc(corr_NoLaser_Spont{v})
    title('Spont - Laser OFF'); caxis([-1 1]); 
    set(gca,'clim',[-1 1],'tickdir','out'); box off; axis square
    subplot(2,3,2); imagesc(corr_Laser_Spont{v})
    title('Spont - Laser ON'); caxis([-1 1]);
    set(gca,'clim',[-1 1],'tickdir','out'); box off; axis square
    subplot(2,3,3); imagesc(corr_Laser_Spont{v} - corr_NoLaser_Spont{v})
    title('Spont: ON - OFF'); caxis([-1 1]);
    set(gca,'clim',[-1 1],'tickdir','out'); box off; axis square
    colorbar
    
    subplot(2,3,4); imagesc(corr_NoLaser_Tone{v})
    title('Tone - Laser OFF'); caxis([-1 1]);
    set(gca,'clim',[-1 1],'tickdir','out'); box off; axis square
    subplot(2,3,5); imagesc(corr_Laser_Tone{v})
    title('Tone - Laser ON'); caxis([-1 1]);
    set(gca,'clim',[-1 1],'tickdir','out'); box off; axis square
    subplot(2,3,6); imagesc(corr_Laser_Tone{v} - corr_NoLaser_Tone{v})
    title('Tone: ON - OFF'); caxis([-1 1]);
    set(gca,'clim',[-1 1],'tickdir','out'); box off; axis square
    colorbar
    set(gcf,'colormap',nanjet);
    set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1200 600]);
set(gcf,'Position',[0 0 1200 600]);
end




%% 
corr_NoLaser_Tone_Dist = [];
corr_Laser_Tone_Dist = [];
corr_NoLaser_Spont_Dist = [];
corr_Laser_Spont_Dist = [];
for v = 1:length(corr_NoLaser_Spont)
    corr_NoLaser_Tone_Dist = [corr_NoLaser_Tone_Dist; corr_NoLaser_Tone{v}(triu(true(size(corr_NoLaser_Tone{v})), 1))];
    corr_Laser_Tone_Dist = [corr_Laser_Tone_Dist; corr_Laser_Tone{v}(triu(true(size(corr_Laser_Tone{v})), 1))];
    corr_NoLaser_Spont_Dist = [corr_NoLaser_Spont_Dist; corr_NoLaser_Spont{v}(triu(true(size(corr_NoLaser_Spont{v})), 1))];
    corr_Laser_Spont_Dist = [corr_Laser_Spont_Dist; corr_Laser_Spont{v}(triu(true(size(corr_Laser_Spont{v})), 1))];

end
figure; 
subplot(2,2,1); histogram(corr_NoLaser_Tone_Dist,25,'EdgeColor','none');
hold on; histogram(corr_Laser_Tone_Dist,25,'FaceColor','k','EdgeColor','none');%alpha(0.5);
line([mean(corr_NoLaser_Tone_Dist) mean(corr_NoLaser_Tone_Dist)],[0 100],'LineStyle','--','Color','b');
line([mean(corr_Laser_Tone_Dist) mean(corr_Laser_Tone_Dist)],[0 100],'LineStyle','--','Color','k');
title('Tone ON (black) vs OFF (blue)'); xlabel('Correlation'); ylabel('# of unit pairs')
box off; axis tight; set(gca,'tickdir','out','xlim',[-1 1],'xtick',[-1 -0.5 0 0.5 1])
subplot(2,2,2); histogram(corr_NoLaser_Spont_Dist,25,'EdgeColor','none');
hold on; histogram(corr_Laser_Spont_Dist,25,'FaceColor','k','EdgeColor','none');%alpha(0.5);
line([nanmean(corr_NoLaser_Spont_Dist) nanmean(corr_NoLaser_Spont_Dist)],[0 100],'LineStyle','--','Color','b');
line([nanmean(corr_Laser_Spont_Dist) nanmean(corr_Laser_Spont_Dist)],[0 100],'LineStyle','--','Color','k');
title('Spont ON (black) vs OFF (blue)');xlabel('Correlation'); ylabel('# of unit pairs')
box off; axis tight; set(gca,'tickdir','out','xlim',[-1 1],'xtick',[-1 -0.5 0 0.5 1])

subplot(2,2,3); histogram(corr_Laser_Tone_Dist - corr_NoLaser_Tone_Dist,25,'EdgeColor','none');
line([mean(corr_Laser_Tone_Dist - corr_NoLaser_Tone_Dist) mean(corr_Laser_Tone_Dist - corr_NoLaser_Tone_Dist)],[0 100],'LineStyle','--','Color','b');
title('Tone Diff (ON - OFF)');xlabel('Correlation difference (ON - OFF)'); ylabel('# of unit pairs')
box off; axis tight; set(gca,'tickdir','out','xlim',[-1 1],'xtick',[-1 -0.5 0 0.5 1])
subplot(2,2,4); histogram(corr_Laser_Spont_Dist - corr_NoLaser_Spont_Dist,25,'EdgeColor','none');
line([nanmean(corr_Laser_Spont_Dist - corr_NoLaser_Spont_Dist) nanmean(corr_Laser_Spont_Dist - corr_NoLaser_Spont_Dist)],[0 100],'LineStyle','--','Color','b');
title('Spont Diff (ON - OFF)');xlabel('Correlation difference (ON - OFF)'); ylabel('# of unit pairs')
box off; axis tight; set(gca,'tickdir','out','xlim',[-1 1],'xtick',[-1 -0.5 0 0.5 1])

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1000 600]);
set(gcf,'Position',[0 0 1000 600]);

cd(FigOut);
print('NoiseCorrelation','-dpdf','-r400')