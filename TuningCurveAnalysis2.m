%% ************************************************************************
%  *****                   1. EXPERIMENT SESSIONS                     *****
%  ************************************************************************
% Load mouse numbers and session numbers for particular set of experimental
% conditions.
opsin = 'chr2';
loc = 'ic';

[MNum, Sesh, TITLE] = loadsessions2(opsin,loc);

FileOutput = 'C:\Users\Jennifer\Documents\MATLAB\TuningCurveAnalysis';
FigOutput = 'C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\ThesisCommittee\Figures\Meeting4';

disp('*******************************************************************');
disp(['1. Experiment sessions loaded: ' TITLE]);
disp('*******************************************************************');
%% ************************************************************************
%  *****                   2. SMOOTH FR                               *****
%  ************************************************************************
% Calculate smoothed firing rate for each cell over highest amplitudes and top preferred frequencies

%Pre-allocate variables
allCELL = cell(1,length(MNum));
FR_LASER = cell(length(MNum),3); 
FR_NOLASER = cell(length(MNum),3);
LOCSoff = cell(1,length(MNum));
LOCSon = cell(1,length(MNum));
CellQ = cell(1,length(MNum));
all_LASERSPIKE2 = cell(1,length(MNum)); 
all_NOLASERSPIKE2 = cell(1,length(MNum));

%Load stimulus parameters
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


for u = 1:length(MNum)  
    %Pull out data files for mice and session numbers of interest
    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    h1 = ls('data\TC003*01.mat');
    idxUSE = find(ismember(str2num(h1(:,13)),Sesh{u}));
    h = h1(idxUSE,:);
    
    %Step 1: Use top 3 amplitudes, but must find the frequencies to use.   
    nfreq = 50; %Number of frequencies to use in smoothes response
    amps = 1;
    [LOCSon{u},LOCSoff{u}] = TC_Select_noGauss(MNum(u),h,nfreq,amps,0);
    
    %Step 2: Separate spikes by tone trial (separated by trial, frequency,
    %amplitude, and laser vs no laser)
    
    Win = 0.24; %Window (in seconds) around start of tone to look for spikes
    
    for n = 1:3

        %Separate out spikes based on if during laser trials vs no laser trials
        LASER_SPIKE = cell(1,size(h,1));
        NOLASER_SPIKE = cell(1,size(h,1));
        for v = 1:size(h,1) 
            match = strfind(h(v,:), '.');
            q = match(end);
            load(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat\TC003_LEFT-0' num2str(n) h(v,6:q-10) '.mat']); 
            CellQ{u}(v) = CellInfo(6);

            SpkTime_NoLaser = SpikeTime(StimOrder_NoLaser,SpikeData,nRep, Win);
            SpkTime_Laser = SpikeTime(StimOrder_Laser,SpikeData,nRep, Win);

            %Select those in bins around top frequencies and top 3 amplitudes
            LASER_SPIKE{v} = SpkTime_Laser(LOCSon{u}(v,:),amps,:);
            NOLASER_SPIKE{v} = SpkTime_NoLaser(LOCSoff{u}(v,:),amps,:);
        end

        %Step 3: Calculate and smooth firing rate for laser and no laser conditions
        SpikeDataLaser = cell(1,size(h,1));        
        SpikeDataNoLaser = cell(1, size(h,1));

        clear LASER_SPIKE2 NOLASER_SPIKE2
        for v = 1:size(h,1)
            for w = 1:numel(LASER_SPIKE{v})
                LASER_SPIKE2{v}{w} = [sort(LASER_SPIKE{v}{w}); w*ones(1,length(LASER_SPIKE{v}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                SpikeDataLaser{v} = [SpikeDataLaser{v} LASER_SPIKE2{v}{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
            end
            FR_LASER{u,n}(v,:) = smoothFRx4(SpikeDataLaser{v},numel(LASER_SPIKE{v})*stimInfo.repeats,0.001,[-Win Win],5);
            %all_LASERSPIKE2{u,n} = LASER_SPIKE2;
            
            for w = 1:numel(NOLASER_SPIKE{v})
                NOLASER_SPIKE2{v}{w} = [sort(NOLASER_SPIKE{v}{w}); w*ones(1,length(NOLASER_SPIKE{v}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                SpikeDataNoLaser{v} = [SpikeDataNoLaser{v} NOLASER_SPIKE2{v}{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
            end
            FR_NOLASER{u,n}(v,:) = smoothFRx4(SpikeDataNoLaser{v},numel(NOLASER_SPIKE{v})*stimInfo.repeats,0.001,[-Win, Win],5);
            %all_NOLASERSPIKE2{u,n} = NOLASER_SPIKE2;
        end
        
    end
        
    allCELL{u} = h;
end

cd(FileOutput);
save(TITLE,'TITLE','MNum','Sesh','allCELL','FR_LASER','FR_NOLASER','CellQ');

disp('*******************************************************************');
disp(['2. FR smoothed for top ' num2str(nfreq) ' frequencies']);
disp('*******************************************************************');

%% ************************************************************************
%  *****               3. FIND TONE RESPONSIVE CELLS                  *****
%  ************************************************************************
% Select good cells and calculate spontaneous firing rate and tone-evoked
% response magnitude for each cell based on first of the TC stimuli

GOODCELL = cell(1,length(MNum));
spontIDX = 185:235; %Indices for spontaneous response
toneIDX = 240:315; %Indices for tone-evoked response
IDX = cell(1,length(MNum));
spontON = cell(1,length(MNum));
spontOFF = cell(1,length(MNum));
aveON = cell(1,length(MNum));
aveOFF = cell(1,length(MNum));
magON = cell(1,length(MNum));
magOFF = cell(1,length(MNum));
idxGOOD = cell(1,length(MNum));

for u = 1:length(MNum)
    %Select good cells
    numcellA = size(allCELL{u},1);
    idxGOOD{u} = NaN(1,numcellA);    
    for j = 1:numcellA
        if CellQ{u}(j) <= 4 && CellQ{u}(j) > 0 && ... %Use only good multi-units and single units
            mean(FR_NOLASER{u,1}(j,240:290)) > (2*std(FR_NOLASER{u,1}(j,145:235)) + mean(FR_NOLASER{u,1}(j,145:235))); %Only use cells with ave tone-evoked FR higher than 2 std's above spontaneous
        
            idxGOOD{u}(j) = 1; %Identify usable cells
            
        end
    end
    IDX{u} = find(~isnan(idxGOOD{u}));
    GOODCELL{u} = allCELL{u}(IDX{u},:);
    
    %Calculate spontaneous firing rate and tone-evoked response magnitude
    numcellB = length(IDX{u});
    spontON{u} = NaN(1,numcellB); spontOFF{u} = NaN(1,numcellB); magON{u} = NaN(1,numcellB);
    magOFF{u} = NaN(1,numcellB); aveON{u} = NaN(1,numcellB); aveOFF{u} = NaN(1,numcellB);
    for j = 1:numcellB   
        spontON{u}(j) = mean(FR_LASER{u,1}(IDX{u}(j),spontIDX)); spontOFF{u}(j) = mean(FR_NOLASER{u,1}(IDX{u}(j),spontIDX));
        aveON{u}(j) = mean(FR_LASER{u,1}(IDX{u}(j),toneIDX)); aveOFF{u}(j) = mean(FR_NOLASER{u,1}(IDX{u}(j),toneIDX));
        magON{u}(j) = aveON{u}(j) - spontON{u}(j); magOFF{u}(j) = aveOFF{u}(j) - spontOFF{u}(j);
    end
    
end

spontOFFall = [spontOFF{:}];
spontONall = [spontON{:}];
magOFFall = [magOFF{:}];
magONall = [magON{:}];


NumberOfCells = length([IDX{:}]);
save(TITLE,'GOODCELL','IDX','magONall','magOFFall','spontONall','spontOFFall','spontIDX','toneIDX','-append')

disp('*******************************************************************');
disp(['3. ' num2str(NumberOfCells) ' cells with tone response']);
disp('*******************************************************************');

%% Make plots
% % ii = 1;
% % % Filters = [0 7 50];
% %  for n = 1%1%1%2:7
% %      cd(['D:\Spikes\M' num2str(MNum(n)) '\'])
% %      close all
% %      numcell = length(IDX{n}); %Number of good cells over both animals
% %      
% % %      for i = 1:length(Filters)
% % %          FR_Laser_M{i} = [FR_LASER{n,i}(ind_cell{n},:)];
% % %          FR_NoLaser_M{i} = [FR_NOLASER{n,i}(ind_cell{n},:)];
% % %      end
% % % 
% % %      %Take mean of laser and no laser across all 3 conditions
% % %      FR_Laser_Mn = [];
% % %      FR_NoLaser_Mn = [];
% % %      for w = 1:length(ind_cell{n})
% % %         FR_Laser_Mn(w,:) = mean([FR_Laser_M{1}(w,:); FR_Laser_M{2}(w,:)]);
% % %         FR_NoLaser_Mn(w,:) = mean([FR_NoLaser_M{1}(w,:); FR_NoLaser_M{2}(w,:); FR_NoLaser_M{3}(w,:)]);
% % %      end
% %     
% %     a = 3; b = 3; %Subplot dimensions (rows vs columns)
% %     nplot = a*b; %number of plots per figure
% %     nfig = ceil(length(IDX{n})./nplot); %Max number of figures needed to plot all cells from mouse n
% % 
% %     for xx = 1:nfig
% %         for j = 1:numcell
% %             if j < (xx*nplot + 1) && j > (xx - 1)*nplot
% %                 figure(xx); 
% %                 subplot(a,b,j - (xx-1)*nplot);
% %                 plot(FR_NOLASER{n,ii}(IDX{n}(j),:),'k')
% %                 hold on; plot(FR_LASER{n,ii}(IDX{n}(j),:),'m')
% %                 MX = max([FR_NOLASER{n,ii}(IDX{n}(j),:) FR_LASER{n,ii}(IDX{n}(j),:)]); %Find max FR out of both laser and no laser
% %                 %area([140 390],[MX MX],'LineStyle','none'); %Location of laser on
% %                 %area([240 290],[MX MX],'LineStyle','none','FaceColor','r'); %Location of tone on
% %                 alpha(0.2)
% %                 hold off;
% %                 title([GOODCELL{1}(j,7:19)]);
% %                 xlabel('Time (ms)'); ylabel('FR (spikes/s)')
% %                 box off
% %                 
% %             end
% %             
% %         end
% %          set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition',[0 0 16 12], 'PaperSize', [16 12], 'color', 'white');
% %          %print('-djpeg','-r400', [num2str(MNum(n)) '_allstim-' num2str(xx)]);
% %      end
% % end
%% ************************************************************************
%  *****              4.CALCULATE AND PLOT SPARSENESS                 *****
%  ************************************************************************
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

%Load tuning curve parameters
load('D:\Code\TuningCurve\TC003_170815LEFT_stimInfo');
fList = stimInfo.index(:,1);
aList = stimInfo.attenuations;
afo5.freqOrder = fList(stimInfo.order)';
afo5.ampOrder = aList(ones(length(stimInfo.order),1))';

nFreq = 50;
amps = 1;
h1 = vertcat(GOODCELL{:});
numcell = size(h1,1);
FR_LASER_ALL = [];
FR_NOLASER_ALL = [];
SidxON = nan(3,numcell);
SidxOFF = nan(3,numcell);


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

        SidxON(n, v) = Sparseness(mFR_ON(n,:),nFreq);  
        SidxOFF(n, v) = Sparseness(mFR_OFF(n,:),nFreq);  
    end

end

descr{1} = {'Laser Onset:'; '-100 ms'};
descr{2} = {'Laser Onset:'; '-20 ms'};
descr{3} = {'Laser Onset:'; '+8 ms'};
figure; 
for n = 1:3
    subplot(3,2,(2*n)-1);
    scatter(SidxOFF(n,:), SidxON(n,:), 25,'filled')
    hold on; line([0 1], [0 1],'Color','k','LineStyle','--')
    hold off; axis square
    xlabel('Sparseness OFF')
    ylabel('Sparseness ON')
    set(gca,'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1])
    [p_sparse(n), h_sparse(n)] = signrank(SidxOFF(n,:), SidxON(n,:));
    text(-0.6,0.5, descr{n});

    subplot(3,2,2*(n));
    bar([nanmean(SidxOFF(n,:)) nanmean(SidxON(n,:))],0.5,'EdgeColor','none');
    hold on; errorbar([nanmean(SidxOFF(n,:)) nanmean(SidxON(n,:))],[nanstd(SidxOFF(n,:))...
        ./sqrt(length(SidxOFF(n,:))) nanstd(SidxON(n,:))./sqrt(length(SidxON(n,:)))],'k','LineStyle','none')
    box off
    ylabel('Sparseness')
    set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out','YLim',[0 1]); axis square
end


set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 1600]);
set(gcf,'Position',[0 0 800 1600]);
suptitle(TITLE)
cd(FigOutput)
print(['sparseness_' TITLE],'-dpdf','-r400')

cd(FileOutput)
save(TITLE,'h_sparse','p_sparse','SidxON','SidxOFF','-append')

disp('*******************************************************************');
disp(['4. Sparseness:         h = ' num2str(h_sparse) ' and p = ' num2str(p_sparse)]);
disp('*******************************************************************');
%% ************************************************************************
%  *****               5. FIND LASER RESPONSIVE CELLS                 *****
%  ************************************************************************
%Identify neurons that show laser effect in spontaneous or tone-evoked

fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

%response for first laser stimulus

%Calculate baseline variability
FRon = cell(1,length(MNum));
FRoff = cell(1,length(MNum));
baseVAR = cell(1,length(MNum));
for u = 1:length(MNum)
    idx = IDX{u};
    FRon{u} = FR_LASER{u,1}(idx,:);
    FRoff{u} = FR_NOLASER{u,1}(idx,:);
    baseVAR{u} = std(FRoff{u}(:,145:235),0,2); %baseline variability
end
LaserThresh = 1; %Threshold
baseVARall = vertcat(baseVAR{:})';

%Calculate spontaneous activity
spontOFFall(isnan(spontOFFall)) = [];
spontONall(isnan(spontONall)) = [];
spontRatio1 = spontONall - (spontOFFall + LaserThresh*baseVARall);
spontRatio2 = spontONall - (spontOFFall - LaserThresh*baseVARall);
SigCellspontUP = find(spontRatio1 > 0); %Decrease activity with laser ON
SigCellspontDOWN = find(spontRatio2 < 0); %Increase activity with laser ON

magOFFall(isnan(magOFFall)) = [];
magONall(isnan(magONall)) = [];
magRatio1 = magONall - (magOFFall + LaserThresh*baseVARall);
magRatio2 = magONall - (magOFFall - LaserThresh*baseVARall);
SigCellmagUP = find(magRatio1 > 0);
SigCellmagDOWN = find(magRatio2 < 0);

BaseFR = []; ToneFR = [];
FR_Laser_M_Norm = []; FR_NoLaser_M_Norm = [];
FR_Laser_M = vertcat(FRon{:});
FR_NoLaser_M = vertcat(FRoff{:});
for j = 1:size(FR_Laser_M,1);
    BaseFR(j) = mean(FR_NoLaser_M(j, 181:230));
    ToneFR(j) = max(FR_NoLaser_M(j,241:290));
    
    FR_Laser_M_Norm(j,:) = (FR_Laser_M(j,:) - BaseFR(j))/(ToneFR(j) - BaseFR(j));
    FR_NoLaser_M_Norm(j,:) = (FR_NoLaser_M(j,:) - BaseFR(j))/(ToneFR(j) - BaseFR(j));
end

%Re-plot ratios to see which have been identified as significant.
figure; 
subplot(3,2,1);
scatter(spontOFFall, spontONall,25,'filled')
hold on; line([0 max([[spontOFFall], [spontONall]])], [0 max([[spontOFFall], [spontONall]])],'Color','k','LineStyle', '--');
%hold on; line([0 31.98], [0 31.98]);
scatter(spontOFFall(SigCellspontDOWN), spontONall(SigCellspontDOWN),25,'m','filled')
scatter(spontOFFall(SigCellspontUP), spontONall(SigCellspontUP),25,'y','filled')
hold off; axis square
set(gca, 'xlim', [0 max([[spontOFFall], [spontONall]])], 'ylim', [0 max([[spontOFFall], [spontONall]])],'TickDir','out')
xlabel('LASER OFF'); ylabel('LASER ON')
title(['Spontaneous activity (Hz)'])
box off

subplot(3,2,2);
scatter(magOFFall, magONall,25,'filled')
hold on; line([0 max([magOFFall, magONall])], [0 max([magOFFall, magONall])],'Color','k','LineStyle', '--');
%hold on; line([0 63.56], [0 63.56]);
scatter(magOFFall(SigCellmagDOWN), magONall(SigCellmagDOWN),25,'m','filled')
scatter(magOFFall(SigCellmagUP), magONall(SigCellmagUP),25,'y','filled')
hold off; axis square
set(gca, 'xlim', [0 max([magOFFall, magONall])], 'ylim', [0 max([magOFFall, magONall])],'TickDir','out')
xlabel('LASER OFF')
title('Tone-evoked response magnitude (Hz)')
box off;

%Plot average time-course for neurons identified as showing significant laser effect
subplot(3,2,3);
plot(mean(FR_NoLaser_M_Norm(SigCellspontDOWN,:),1),'k','LineWidth',2)
hold on; plot(mean(FR_Laser_M_Norm(SigCellspontDOWN,:),1),'m','LineWidth',2); hold off;
ylabel('Normalized FR'); title(['N = ' num2str(length(SigCellspontDOWN))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:100:480,'XTickLabel',[-0.1:0.1:0.3], ...
    'YLim', [-0.2 1.2],'YTick', 0:0.2:1)
subplot(3,2,5);
plot(mean(FR_NoLaser_M_Norm(SigCellspontUP,:),1),'k','LineWidth',2)
hold on; plot(mean(FR_Laser_M_Norm(SigCellspontUP,:),1),'y','LineWidth',2); hold off;
ylabel('Normalized FR'); xlabel('Time (s)'); title(['N = ' num2str(length(SigCellspontUP))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:100:480,'XTickLabel',[-0.1:0.1:0.3], ...
    'YLim', [-0.2 1.2],'YTick', 0:0.2:1)

subplot(3,2,4);
plot(mean(FR_NoLaser_M_Norm(SigCellmagDOWN,:),1),'k','LineWidth',2)
hold on; plot(mean(FR_Laser_M_Norm(SigCellmagDOWN,:),1),'m','LineWidth',2); hold off;
title(['N = ' num2str(length(SigCellmagDOWN))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:100:480,'XTickLabel',[-0.1:0.1:0.3], ...
    'YLim', [-0.2 1.2],'YTick', 0:0.2:1)
subplot(3,2,6);  
plot(mean(FR_NoLaser_M_Norm(SigCellmagUP,:),1),'k','LineWidth',2)
hold on; plot(mean(FR_Laser_M_Norm(SigCellmagUP,:),1),'y','LineWidth',2); hold off;
xlabel('Time (s)'); title(['N = ' num2str(length(SigCellmagUP))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:100:480,'XTickLabel',[-0.1:0.1:0.3], ...
    'YLim', [-0.2 1.2],'YTick', 0:0.2:1)

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1000 1200]);
set(gcf,'Position',[0 0 1000 1200]);
suptitle(TITLE)
cd(FigOutput)
print(['laser_' TITLE],'-dpdf','-r400')

%Identify which cells show a laser effect in either direction
IDXall = horzcat(IDX{:});
GOODCELLall = vertcat(GOODCELL{:});
IDXup = IDXall(unique([SigCellmagUP SigCellspontUP]));
GOODCELLup = GOODCELLall(unique([SigCellmagUP SigCellspontUP]),:);
IDXdown = IDXall(unique([SigCellmagDOWN SigCellspontDOWN]));
GOODCELLdown = GOODCELLall(unique([SigCellmagDOWN SigCellspontDOWN]),:);
IDXsig = IDXall(unique([SigCellmagUP SigCellspontUP SigCellmagDOWN SigCellspontDOWN]));
GOODCELLsig = GOODCELLall(unique([SigCellmagUP SigCellspontUP SigCellmagDOWN SigCellspontDOWN]),:);

cd(FileOutput)
save(TITLE,'SigCellmagUP','SigCellmagDOWN','SigCellspontUP','SigCellspontDOWN','-append')

disp('*******************************************************************');
disp(['5. ' num2str(length(IDXsig)) ' show change'])
disp('*******************************************************************');

%% ************************************************************************
%  *****              6.CALCULATE AND PLOT SPARSENESS                 *****
%  ************************************************************************
%Calculate sparseness for cells that show change due to laser
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

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
h1 = GOODCELLsig;
numcell = size(h1,1);
FR_LASER_ALL = [];
FR_NOLASER_ALL = [];
SidxON = nan(3,numcell);
SidxOFF = nan(3,numcell);


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

    subplot(3,2,2*(n));
    bar([nanmean(LASERSidxOFF(n,:)) nanmean(LASERSidxON(n,:))],0.5,'EdgeColor','none');
    hold on; errorbar([nanmean(LASERSidxOFF(n,:)) nanmean(LASERSidxON(n,:))],[nanstd(LASERSidxOFF(n,:))...
        ./sqrt(length(LASERSidxOFF(n,:))) nanstd(LASERSidxON(n,:))./sqrt(length(LASERSidxON(n,:)))],'k','LineStyle','none','LineWidth',2)
    box off
    ylabel('Sparseness')
    set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out','YLim',[0 1]); axis square
end


set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 1600]);
set(gcf,'Position',[0 0 800 1600]);
suptitle([TITLE ': Affected cells'])
cd(FigOutput)
print(['LASERsparseness_' TITLE],'-dpdf','-r400')

cd(FileOutput)
save(TITLE,'h_LASERsparse','p_LASERsparse','LASERSidxON','LASERSidxOFF','-append')

disp('*******************************************************************');
disp(['6. Sparseness of affected cells: h = ' num2str(h_sparse) ' and p = ' num2str(p_sparse)]);
disp('*******************************************************************');

%% ************************************************************************
%  *****                         7.LINEAR FITS                        *****
%  ************************************************************************
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

cd('D:\Spikes')
onsets = {'-100 ms','-20 ms', '+8 ms'};
use = GOODCELLall(SigCellmagDOWN,:);
fitted = cell(1,3);
linearparams{1} = NaN(size(use,1),3);
figure;
for n = 1:3
    for i = 1:size(use,1)
        q = find(use(i,:) == '.');
        load([use(i,1:q - 2) num2str(n) '.mat'])
        %Normalize firing rate
        Ton = sort(TCon.TCmat{1});
        Toff = sort(TCoff.TCmat{1});
        ToffNorm = (Toff - min(Toff))./(max(Toff) - min(Toff));
        TonNorm = (Ton - min(Toff))./(max(Toff) - min(Toff));

        %Calculate linear fits
        mdl = fitlm(ToffNorm,TonNorm);
        linearparams{1}(i,n) = mdl.Coefficients{2,1};
        linearparams{1}(i,n+3) = mdl.Coefficients{1,1};
        fitted{n}(i,:) = mdl.Coefficients{2,1}*[0:0.05:1]+mdl.Coefficients{1,1};
    end
    
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
LinearFits.magDOWN = fitted;
suptitle('Feedback Affected Cells: Response Mag DECREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOutput)
print('LinFit_magDOWN','-dpdf','-r400')


cd('D:\Spikes')
use = GOODCELLall(SigCellmagUP,:);
fitted = cell(1,3);
linearparams{2} = NaN(size(use,1),3);
figure;
for n = 1:3
    for i = 1:size(use,1)
        q = find(use(i,:) == '.');
        load([use(i,1:q - 2) num2str(n) '.mat'])
        %Normalize firing rate
        Ton = sort(TCon.TCmat{1});
        Toff = sort(TCoff.TCmat{1});
        ToffNorm = (Toff - min(Toff))./(max(Toff) - min(Toff));
        TonNorm = (Ton - min(Toff))./(max(Toff) - min(Toff));

        %Calculate linear fits
        mdl = fitlm(ToffNorm,TonNorm);
        linearparams{2}(i,n) = mdl.Coefficients{2,1};
        linearparams{2}(i,n+3) = mdl.Coefficients{1,1};
        fitted{n}(i,:) = mdl.Coefficients{2,1}*[0:0.05:1]+mdl.Coefficients{1,1};
    end
    
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
LinearFits.magUP = fitted;
suptitle('Feedback Affected Cells: Response Mag INCREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOutput)
print('LinFit_magUP','-dpdf','-r400')

cd('D:\Spikes')
use = GOODCELLall(SigCellspontDOWN,:);
linearparams{3} = NaN(size(use,1),3);
fitted = cell(1,3);
figure;
for n = 1:3
    for i = 1:size(use,1)
        q = find(use(i,:) == '.');
        load([use(i,1:q - 2) num2str(n) '.mat'])
        %Normalize firing rate
        Ton = sort(TCon.TCmat{1});
        Toff = sort(TCoff.TCmat{1});
        ToffNorm = (Toff - min(Toff))./(max(Toff) - min(Toff));
        TonNorm = (Ton - min(Toff))./(max(Toff) - min(Toff));
        
        %Calculate linear fits
        mdl = fitlm(ToffNorm,TonNorm);
        linearparams{3}(i,n) = mdl.Coefficients{2,1};
        linearparams{3}(i,n+3) = mdl.Coefficients{1,1};
        fitted{n}(i,:) = mdl.Coefficients{2,1}*[0:0.05:1]+mdl.Coefficients{1,1};
    end
    
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
LinearFits.spontDOWN = fitted;
suptitle('Feedback Affected Cells: Response Spont DECREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOutput)
print('LinFit_spontDOWN','-dpdf','-r400')

cd('D:\Spikes')
use = GOODCELLall(SigCellspontUP,:);
linearparams{4} = NaN(size(use,1),3);
fitted = cell(1,3);
figure;
for n = 1:3
    for i = 1:size(use,1)
        q = find(use(i,:) == '.');
        load([use(i,1:q - 2) num2str(n) '.mat'])
        %Normalize firing rate
        Ton = sort(TCon.TCmat{1});
        Toff = sort(TCoff.TCmat{1});
        ToffNorm = (Toff - min(Toff))./(max(Toff) - min(Toff));
        TonNorm = (Ton - min(Toff))./(max(Toff) - min(Toff));
        
        %Calculate linear fits
        mdl = fitlm(ToffNorm,TonNorm);
        linearparams{4}(i,n) = mdl.Coefficients{2,1};
        linearparams{4}(i,n+3) = mdl.Coefficients{1,1};
        fitted{n}(i,:) = mdl.Coefficients{2,1}*[0:0.05:1]+mdl.Coefficients{1,1};
    end
    
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
LinearFits.spontUP = fitted;

suptitle('Feedback Affected Cells: Response Spont INCREASE')
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 500]);
set(gcf,'Position',[0 0 1600 500]);
cd(FigOutput)
print('LinFit_spontUP','-dpdf','-r400')

cd(FileOutput)
save(TITLE,'LinearFits','linearparams','-append')

%% ************************************************************************
%  *****                        7.STRF ANALYSIS                       *****
%  ************************************************************************

%parameters
its = 1; %iterations for random STA

sigma_b = 1.5;
p = .05;	%p-value
tC = 2.3; 	%threshold (for Cluster test)
LaserOn1 = 0.5;
LaserDur = 0.25;
stimdur = 300;
win = [0 0.25; 0.25 0.5; 0.5 0.75; 0.75 1];
STAwin = 0.1;

%Load stimulus parameters
load('DRC001-01','params');
randparams = params;
times = [-0.1:0.005:0];
freqs = params.freqs;
nFiles = 8;

temp = vertcat(allCELL{:});
Qcell = temp(horzcat(CellQ{:}) > 0 & horzcat(CellQ{:}) < 5,:);
strfFR = NaN(size(Qcell,1),3);
for u = 20%1:size(Qcell,1)
    q = find(Qcell(u,:) == '.');
    load(['DRC001-' Qcell(u,7:q-10) '.mat'])
    load(['DRC001_LEFT-01-' Qcell(u,7:q-10) '.mat'])
    
%     %%%Load spike file here (also append cellInfo with STRFclust)%%%
    spikes = SpikeData(3,:);
    spikes(spikes >= stimdur) = [];

    randSTRFall = zeros(size(STAon,1),size(STAon,2),its);
    fprintf('Generating random STA(s)...\n')
        for i = 1:its
            fprintf('Iteration %02d/%02d...\n',i,its);
            % make a stim shuffle index
            idx = randperm(size(params.dbs,2));
            randparams.dbs = params.dbs(:,idx);
            % make a new sta from shuffled stim
            randSTRFall(:,:,i) = genSTRF(spikes,randparams,STAwin);
        end

        % get the mean matrix of the permutation - this is the noise:
        randSTRF = mean(randSTRFall,3);  
        

        STRFclustON.POS = calcSTRFcluster(STAon,randSTRF,sigma_b,p,tC);
        STRFclustON.NEG = calcSTRFcluster(-STAon,-randSTRF,sigma_b,p,tC);
        STRFclustOFF.POS = calcSTRFcluster(STAoff3,randSTRF,sigma_b,p,tC);
        STRFclustOFF.NEG = calcSTRFcluster(-STAoff3,-randSTRF,sigma_b,p,tC);
        
        %Plot masks for positive lobe just to check if cluster test seems
        %reasonable
        MM = max(max([STAon;STAoff3]));
        mm = min(min([STAon;STAoff3]));
        subplot(2,2,1); plotSTA([-0.1:0.005:0],params.freqs/1000,STAon,1,[mm MM]);

        subplot(2,2,3); plotSTA([-0.1:0.005:0],params.freqs/1000,STAoff3,1,[mm MM]);
        subplot(2,2,4); imagesc(STRFclustOFF.POS.tCi.ClusTest)
        set(gca,'YDir','normal')
        subplot(2,2,2); imagesc(STRFclustON.POS.tCi.ClusTest)
        set(gca,'YDir','normal')
        %pause(0.1)
        
      
        
        if STRFclustOFF.POS.info.Clusterdata(end) == 1 && STRFclustON.POS.info.Clusterdata(end) == 1 
            [ClusMask, STRFclustON.POS.params.OFFmaskmean] = matchSTRFclust(STRFclustOFF.POS.tCi.ClusMask,STRFclustON.POS.tCi.ClusMask,STAon);
            STRFclustON.POS.tCi.ClusMaskOrig = STRFclustON.POS.tCi.ClusMask;
            STRFclustON.POS.tCi.ClusMask = ClusMask;
        end
                           
        
        if STRFclustOFF.NEG.info.Clusterdata(end) == 1 && STRFclustON.NEG.info.Clusterdata(end) == 1 
            [ClusMask, STRFclustON.NEG.params.OFFmaskmean] = matchSTRFclust(STRFclustOFF.NEG.tCi.ClusMask,STRFclustON.NEG.tCi.ClusMask,STAon);
            STRFclustON.NEG.tCi.ClusMaskOrig = STRFclustON.NEG.tCi.ClusMask;
            STRFclustON.NEG.tCi.ClusMask = ClusMask;
        end
                
        STRFclustOFF.POS = calcSTRFparams(STRFclustOFF.POS,STAoff3,times,freqs);
        STRFclustON.POS = calcSTRFparams(STRFclustON.POS,STAon,times,freqs);
        STRFclustOFF.NEG = calcSTRFparams(STRFclustOFF.NEG,STAoff3,times,freqs);
        STRFclustON.NEG = calcSTRFparams(STRFclustON.NEG,STAon,times,freqs);  
        
        %save(['DRC001-' Qcell(u,7:q-10) '.mat'],'STRFclustON','STRFclustOFF','-append')
        
        for ii = 1:nFiles

            load(['SpikeMat\DRC001_LEFT-0' num2str(ii) '-' Qcell(u,7:q-10) '.mat']); %Load spike times 

            spikes = SpikeData(3,:);
            spikes(spikes >= stimdur) = [];
            SpikeTmod = mod(spikes,1);
            SpikeFR(ii,:) = smoothFRx4(SpikeTmod,stimdur,0.001,[0 1],5);
        end
        tempFR = mean(SpikeFR,1);
        strfFR(u,1) = mean(tempFR(500:750)); %Laser ON firing rate
        strfFR(u,2) = mean(tempFR(1:499)); %Laser OFF firing rate
        strfFR(u,3) = std(tempFR(1:499));     

end

delta = strfFR(:,1) - (strfFR(:,2) - 2*strfFR(:,3));
strfUP = find(delta > 0);
strfcellUP = Qcell(strfUP,:);
delta = (strfFR(:,2) - 2*strfFR(:,3)) - strfFR(:,1);
strfDOWN = find(delta > 0);
strfcellDOWN = Qcell(strfDOWN,:);

%cd(FileOutput)
%save(TITLE,'strfcellDOWN','strfcellUP','-append')

