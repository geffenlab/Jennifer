%% ************************************************************************
%  *****                   1. EXPERIMENT SESSIONS                     *****
%  ************************************************************************
% Load mouse numbers and session numbers for particular set of experimental
% conditions.
opsin = 'archt';
mouseline = 'camk2';
loc = 'ic';
cond = 'awake';

[MNum, Sesh, TITLE] = loadsessions(opsin,mouseline,loc,cond);

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
FR_LASER = cell(1,length(MNum)); 
FR_NOLASER = cell(1,length(MNum));
LOCSoff = cell(1,length(MNum));
LOCSon = cell(1,length(MNum));
CellQ = cell(1,length(MNum));
all_LASERSPIKE2 = cell(1,length(MNum)); 
all_NOLASERSPIKE2 = cell(1,length(MNum));

afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
Win = 0.24; %Window (in seconds) around start of tone to look for spikes
trialdur = afo5.interval+afo5.soundLen;
dur = trialdur*length(afo5.freqOrder);
Time = 0:trialdur:dur; %Each repetition is 400 seconds long, each trial is 500ms long
Time = Time(1:end-1);

StimOrder_Laser = [Time(2:2:end); afo5.freqOrder(2:2:end); afo5.ampOrder(2:2:end)];
StimOrder_NoLaser = [Time(1:2:end); afo5.freqOrder(1:2:end); afo5.ampOrder(1:2:end)];


for u = 1:length(MNum)  
    %Pull out data files for mice and session numbers of interest
    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    h1 = ls('data\TC3-1*_laser.mat');
    idxUSE = find(ismember(str2num(h1(:,15)),Sesh{u}));
    h = h1(idxUSE,:);
    
    
    %Step 1: Use top 3 amplitudes, but must find the frequencies to use.   
    nfreq = 7; %Number of frequencies to use in smoothes response
    amps = [6:8];
    [LOCSon{u},LOCSoff{u}] = TC_Select_noGauss(MNum(u),h,nfreq,amps,0);
    
    %Step 2: Separate spikes by tone trial (separated by trial, frequency,
    %amplitude, and laser vs no laser)


    %Separate out spikes based on if during laser trials vs no laser trials
    LASER_SPIKE = cell(1,size(h,1));
    NOLASER_SPIKE = cell(1,size(h,1));
    for v = 1:size(h,1) 
        match = strfind(h(v,:), '0');
        q = match(end);
        load(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat\R407F'  h(v,5:q) '.mat']); 
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
            FR_LASER{u}(v,:) = smoothFRx4(SpikeDataLaser{v},numel(LASER_SPIKE{v}),0.001,[-Win Win],5);
            all_LASERSPIKE2{u} = LASER_SPIKE2;
            
            for w = 1:numel(NOLASER_SPIKE{v})
                NOLASER_SPIKE2{v}{w} = [sort(NOLASER_SPIKE{v}{w}); w*ones(1,length(NOLASER_SPIKE{v}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                SpikeDataNoLaser{v} = [SpikeDataNoLaser{v} NOLASER_SPIKE2{v}{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
            end
            FR_NOLASER{u}(v,:) = smoothFRx4(SpikeDataNoLaser{v},numel(NOLASER_SPIKE{v}),0.001,[-Win, Win],5);
            all_NOLASERSPIKE2{u} = NOLASER_SPIKE2;
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
% Select good cells and calculate spontaneous firing rate and tone-evoked response magnitude for each cell

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
            mean(FR_NOLASER{u}(j,240:290)) > (2*std(FR_NOLASER{u}(j,145:235)) + mean(FR_NOLASER{u}(j,145:235))); %Only use cells with ave tone-evoked FR higher than 2 std's above spontaneous
        
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
        spontON{u}(j) = mean(FR_LASER{u}(IDX{u}(j),spontIDX)); spontOFF{u}(j) = mean(FR_NOLASER{u}(IDX{u}(j),spontIDX));
        aveON{u}(j) = mean(FR_LASER{u}(IDX{u}(j),toneIDX)); aveOFF{u}(j) = mean(FR_NOLASER{u}(IDX{u}(j),toneIDX));
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

%% ************************************************************************
%  *****                  4. FIND USABLE CELLS                        *****
%  ************************************************************************
% Categorize each neuron based on spontaneous activity (to remove units
% with ONLY laser onset/offest responses) and remove cells with laser
% artifact.

%Identify neurons that show laser effect in each of 3 categories from
%previous section
FRon = cell(1,length(MNum));
FRoff = cell(1,length(MNum));
baseVAR = cell(1,length(MNum));
OnsetResp1 = cell(1,length(MNum));
OnsetResp2 = cell(1,length(MNum));
OffsetResp = cell(1,length(MNum));
for n = 1:length(MNum)
    if ~isempty(IDX{n})
        FRon{n} = FR_LASER{n}(IDX{n},:);
        FRoff{n} = FR_NOLASER{n}(IDX{n},:);
        baseVAR{n} = std(FRoff{n}(:,spontIDX),0,2); %baseline variability
        OnsetResp1{n} = mean(FRon{n}(:,140:180),2);
        OnsetResp2{n} = mean(FRon{n}(:,200:240),2);
        OffsetResp{n} = mean(FRon{n}(:,390:430),2);
    end
end
Thresh = 3; %Threshold
baseVARall = vertcat(baseVAR{:})';
OnsetResp1all = vertcat(OnsetResp1{:})';
OnsetResp2all = vertcat(OnsetResp2{:})';
OffsetRespall = vertcat(OffsetResp{:})';
nCell = length(baseVARall);
GOODCELLall = strvcat(GOODCELL{:});

%Categorize cells
CellCat = NaN(1,nCell);
for i = 1:nCell
    if OnsetResp2all(i) > spontOFFall(i) + Thresh*baseVARall(i) || OnsetResp2all(i) < spontOFFall(i) - Thresh*baseVARall(i) %Must have significant difference 25 ms before tone onset
        CellCat(i) = 1; % SUSTAINED LASER RESPONSE
    elseif OnsetResp1all(i) > spontOFFall(i) + Thresh*baseVARall(i) || OnsetResp1all(i) < spontOFFall(i) - Thresh*baseVARall(i)
        CellCat(i) = 2; % ONSET RESPONSE 
    elseif OffsetRespall(i) > spontOFFall(i) + Thresh*baseVARall(i) || OffsetRespall(i) < spontOFFall(i) - Thresh*baseVARall(i)
        CellCat(i) = 3; % OFFSET ONLY RESPONSE 
    elseif spontONall(i) < spontOFFall(i) + Thresh*baseVARall(i) && spontONall(i) > spontOFFall(i) - Thresh*baseVARall(i)
        CellCat(i) = 0; % NO LASER RESPONSE 
    end    
end

%Select only units that are category 0 and category 1
badIDX = find(CellCat == 2 | CellCat == 3 | isnan(CellCat));
GOODCELLall(badIDX,:) = [];

spontOFFall2 = spontOFFall;
spontOFFall2(badIDX) = [];
spontONall2 = spontONall;
spontONall2(badIDX) = [];

magOFFall2 = magOFFall;
magOFFall2(badIDX) = [];
magONall2 = magONall;
magONall2(badIDX) = [];

NumberOfCellsUse = length(spontOFFall2);

save(TITLE,'GOODCELLall','magONall2','magOFFall2','spontONall2','spontOFFall2','Thresh','-append')

disp('*******************************************************************');
disp(['4. ' num2str(NumberOfCellsUse) ' usable cells']);
disp('*******************************************************************');
%% ************************************************************************
%  *****          5. PLOT TONE-EVOKED RESPONSE MAGNITUDE              *****
%  ************************************************************************
% Use category 0 and category 1 cells to plot tone-evoked response magnitude and spontaneous activity

% Need to adjust limits if combining anesthetized and awake onto a single
% figure.

fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

%figure('rend','painters','pos',[10 400 2000 400]); 
figure;
subplot(1,4,1); scatter(spontOFFall2, spontONall2,25,'filled')
LIM = max([spontOFFall2, spontONall2]);
LIM = 50;
hold on; line([0 LIM], [0 LIM],'Color','k','LineStyle','--');
hold off;
set(gca,'TickDir','out'); axis square
set(gca,'xlim',[0 LIM],'ylim',[0 LIM])
xlabel('LASER OFF')
ylabel('LASER ON')
title('Spontaneous activity (Hz)')

subplot(1,4,2); bar([mean(spontOFFall2) mean(spontONall2)],0.5,'EdgeColor','none');
hold on; errorbar([mean(spontOFFall2) mean(spontONall2)],[std(spontOFFall2)./sqrt(length(spontOFFall2)) std(spontONall2)./sqrt(length(spontONall2))],'k','LineStyle','none','LineWidth',2)
box off
ylabel('Firing rate (Hz)')
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
set(gca,'YLim',[0 15])
[p_spont, h_spont] = signrank(spontOFFall2, spontONall2);

subplot(1,4,3); scatter(magOFFall2, magONall2,25,'filled')
LIM = max([magOFFall2, magONall2]);
LIM = 50;
hold on; line([0 LIM], [0 LIM], 'Color','k','LineStyle','--');
hold off;
set(gca,'TickDir','out'); axis square
set(gca,'xlim',[0 LIM],'ylim',[0 LIM])
xlabel('LASER OFF')
ylabel('LASER ON')
title('Tone-evoked response magnitude (Hz)')

subplot(1,4,4); bar([mean(magOFFall2) mean(magONall2)],0.5,'EdgeColor','none');
hold on; errorbar([mean(magOFFall2) mean(magONall2)],[std(magOFFall2)./sqrt(length(magOFFall2)) std(magONall2)./sqrt(length(magONall2))],'k','LineStyle','none','LineWidth',2)
box off
ylabel('Firing rate (Hz)')
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
set(gca,'YLim',[0 30])
[p_mag, h_mag] = signrank(magOFFall2, magONall2);

suptitle(TITLE);
cd(FileOutput);
save(TITLE,'h_mag','p_mag','h_spont','p_spont','-append')

cd(FigOutput)
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 400]);
set(gcf,'Position',[0 0 1600 400]);
print(TITLE,'-dpdf','-r400')


disp('*******************************************************************');
disp(['5. Response magnitude: h = ' num2str(h_mag) ' and p = ' num2str(p_mag)]);
disp(['   Spontaneous:        h = ' num2str(h_spont) ' and p = ' num2str(p_spont)]);
disp('*******************************************************************');
%% ************************************************************************
%  *****              6.CALCULATE AND PLOT SPARSENESS                 *****
%  ************************************************************************
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

nFreq = 50;
numcell = size(GOODCELLall,1);
FR_LASER_ALL = [];
FR_NOLASER_ALL = [];

h1 = GOODCELLall;
afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
Win = 0.24; %Window (in seconds) around start of tone to look for spikes
trialdur = afo5.interval+afo5.soundLen;
dur = trialdur*length(afo5.freqOrder);
Time = 0:trialdur:dur; %Each repetition is 400 seconds long, each trial is 500ms long
Time = Time(1:end-1);

StimOrder_Laser = [Time(2:2:end); afo5.freqOrder(2:2:end); afo5.ampOrder(2:2:end)];
StimOrder_NoLaser = [Time(1:2:end); afo5.freqOrder(1:2:end); afo5.ampOrder(1:2:end)];


SpkTime_LaserAll = cell(size(h1,1),nFreq);
SpkTime_NoLaserAll = cell(size(h1,1),nFreq);

for v = 1:size(h1,1) %for filter 3, v = 47 peak is at index 2

    match = strfind(h1(v,:), '_');
    q = match(end);
    load(['D:\Spikes\M' h1(v,10:13) '\SpikeMat\R407F'  h1(v,5:q-1) '.mat']); 

    SpkTime_NoLaser = SpikeTime(StimOrder_NoLaser,SpikeData,nRep, Win);
    SpkTime_Laser = SpikeTime(StimOrder_Laser,SpikeData,nRep, Win);
    for j = 1:nFreq
        SpkTime_LaserAll{v,j} = SpkTime_Laser(j,6:8,:);
        SpkTime_NoLaserAll{v,j} = SpkTime_NoLaser(j,6:8,:);
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
        FR_LASER_ALL(j,v,:) = smoothFRx4(SpikeDataLaserAll{v,j},numel(LASER_SPIKE_ALL{v,j}),0.001,[-Win Win],5);

        %For no laser trials
        for w = 1:numel(SpkTime_NoLaserAll{v})
            NOLASER_SPIKE_ALL{v,j}{w} = [sort(SpkTime_NoLaserAll{v,j}{w}); w*ones(1,length(SpkTime_NoLaserAll{v,j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
            SpikeDataNoLaserAll{v,j} = [SpikeDataNoLaserAll{v,j} NOLASER_SPIKE_ALL{v,j}{w}];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
        end
        FR_NOLASER_ALL(j,v,:) = smoothFRx4(SpikeDataNoLaserAll{v,j},numel(NOLASER_SPIKE_ALL{v,j}),0.001,[-Win Win],5);

    end
end
   

%Calculate Sparseness index for each cell
SidxON = nan(1,numcell);
SidxOFF = nan(1,numcell);
for v = 1:numcell
   for j = 1:nFreq
       mFR_ON(j) = mean(FR_LASER_ALL(j,v,241:290));
       mFR_OFF(j) = mean(FR_NOLASER_ALL(j,v,241:290));
   end
       
    SidxON(v) = Sparseness(mFR_ON,nFreq);  
    SidxOFF(v) = Sparseness(mFR_OFF,nFreq);  
end

figure; subplot(1,2,1);
scatter(SidxOFF, SidxON, 25,'filled')
hold on; line([0 1], [0 1],'Color','k','LineStyle','--')
hold off; axis square
xlabel('Sparseness OFF')
ylabel('Sparseness ON')
title('Sparseness')
set(gca,'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1])
[p_sparse, h_sparse] = signrank(SidxOFF, SidxON);

subplot(1,2,2);
bar([nanmean(SidxOFF) nanmean(SidxON)],0.5,'EdgeColor','none');
hold on; errorbar([nanmean(SidxOFF) nanmean(SidxON)],[nanstd(SidxOFF)./sqrt(length(SidxOFF)) nanstd(SidxON)./sqrt(length(SidxON))],'k','LineStyle','none','LineWidth',2)
box off
ylabel('Sparseness')
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out','YLim',[0 1]); axis square

suptitle(TITLE)
cd(FigOutput)
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 400]);
set(gcf,'Position',[0 0 800 400]);
print(['sparseness_' TITLE],'-dpdf','-r400')

cd(FileOutput)
save(TITLE,'h_sparse','p_sparse','SidxON','SidxOFF','-append')

disp('*******************************************************************');
disp(['6. Sparseness:         h = ' num2str(h_sparse) ' and p = ' num2str(p_sparse)]);
disp('*******************************************************************');

%% ************************************************************************
%  *****               7. FIND LASER RESPONSIVE CELLS                 *****
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
    a = ismember(GOODCELL{u},GOODCELLall,'rows');
    idx = IDX{u}(a);
    if ~isempty(idx)
    FRon{u} = FR_LASER{u}(idx,:);
    FRoff{u} = FR_NOLASER{u}(idx,:);
    baseVAR{u} = std(FRoff{u}(:,145:235),0,2); %baseline variability
    end
end
LaserThresh = 1; %Threshold
baseVARall = vertcat(baseVAR{:})';

%Calculate spontaneous activity
spontRatio1 = spontONall2 - (spontOFFall2 + LaserThresh*baseVARall);
spontRatio2 = spontONall2 - (spontOFFall2 - LaserThresh*baseVARall);
SigCellspontUP = find(spontRatio1 > 0); %Decrease activity with laser ON
SigCellspontDOWN = find(spontRatio2 < 0); %Increase activity with laser ON

magRatio1 = magONall2 - (magOFFall2 + LaserThresh*baseVARall);
magRatio2 = magONall2 - (magOFFall2 - LaserThresh*baseVARall);
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
scatter(spontOFFall2, spontONall2,25,'filled')
hold on; line([0 max([[spontOFFall2], [spontONall2]])], [0 max([[spontOFFall2], [spontONall2]])],'Color','k','LineStyle', '--');
%hold on; line([0 31.98], [0 31.98]);
scatter(spontOFFall2(SigCellspontDOWN), spontONall2(SigCellspontDOWN),25,'m','filled')
scatter(spontOFFall2(SigCellspontUP), spontONall2(SigCellspontUP),25,'y','filled')
hold off; axis square
set(gca, 'xlim', [0 max([[spontOFFall2], [spontONall2]])], 'ylim', [0 max([[spontOFFall2], [spontONall2]])],'TickDir','out')
xlabel('LASER OFF'); ylabel('LASER ON')
title(['Spontaneous activity (Hz)'])
box off

subplot(3,2,2);
scatter(magOFFall2, magONall2,25,'filled')
hold on; line([0 max([magOFFall2, magONall2])], [0 max([magOFFall2, magONall2])],'Color','k','LineStyle', '--');
%hold on; line([0 63.56], [0 63.56]);
scatter(magOFFall2(SigCellmagDOWN), magONall2(SigCellmagDOWN),25,'m','filled')
scatter(magOFFall2(SigCellmagUP), magONall2(SigCellmagUP),25,'y','filled')
hold off; axis square
set(gca, 'xlim', [0 max([magOFFall2, magONall2])], 'ylim', [0 max([magOFFall2, magONall2])],'TickDir','out')
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
%print(['laser_' TITLE],'-dpdf','-r400')

%Identify which cells show a laser effect in either direction
IDXall = horzcat(IDX{:});
GOODCELLall = strvcat(GOODCELL{:});
IDXup = IDXall(unique([SigCellmagUP SigCellspontUP]));
GOODCELLup = GOODCELLall(unique([SigCellmagUP SigCellspontUP]),:);
IDXdown = IDXall(unique([SigCellmagDOWN SigCellspontDOWN]));
GOODCELLdown = GOODCELLall(unique([SigCellmagDOWN SigCellspontDOWN]),:);
IDXsig = IDXall(unique([SigCellmagUP SigCellspontUP SigCellmagDOWN SigCellspontDOWN]));
GOODCELLsig = GOODCELLall(unique([SigCellmagUP SigCellspontUP SigCellmagDOWN SigCellspontDOWN]),:);

disp('*******************************************************************');
disp(['7. ' num2str(length(IDXsig)) ' show change'])
disp('*******************************************************************');

%% ************************************************************************
%  *****              8.CALCULATE AND PLOT SPARSENESS                 *****
%  ************************************************************************
%Calculate sparseness for cells that show change due to laser
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

%Load tuning curve parameters
afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters

nFreq = 50;
amps = 1;
h1 = GOODCELLsig;
numcell = size(h1,1);
FR_LASER_ALL = [];
FR_NOLASER_ALL = [];
LASERSidxON = nan(1,numcell);
LASERSidxOFF = nan(1,numcell);


afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
Win = 0.24; %Window (in seconds) around start of tone to look for spikes
trialdur = afo5.interval+afo5.soundLen;
dur = trialdur*length(afo5.freqOrder);
Time = 0:trialdur:dur; %Each repetition is 400 seconds long, each trial is 500ms long
Time = Time(1:end-1);

StimOrder_Laser = [Time(2:2:end); afo5.freqOrder(2:2:end); afo5.ampOrder(2:2:end)];
StimOrder_NoLaser = [Time(1:2:end); afo5.freqOrder(1:2:end); afo5.ampOrder(1:2:end)];


SpkTime_LaserAll = cell(size(h1,1),nFreq);
SpkTime_NoLaserAll = cell(size(h1,1),nFreq);

for v = 1:size(h1,1) 

    match = strfind(h1(v,:), '_');
    q = match(end);
    load(['D:\Spikes\M' h1(v,10:13) '\SpikeMat\R407F'  h1(v,5:q-1) '.mat']); 

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
        FR_LASER_ALL(j,v,:) = smoothFRx4(SpikeDataLaserAll{v,j},numel(LASER_SPIKE_ALL{v,j}),0.001,[-Win Win],5);

        %For no laser trials
        for w = 1:numel(SpkTime_NoLaserAll{v})
            NOLASER_SPIKE_ALL{v,j}{w} = [sort(SpkTime_NoLaserAll{v,j}{w}); w*ones(1,length(SpkTime_NoLaserAll{v,j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
            SpikeDataNoLaserAll{v,j} = [SpikeDataNoLaserAll{v,j} NOLASER_SPIKE_ALL{v,j}{w}];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
        end
        FR_NOLASER_ALL(j,v,:) = smoothFRx4(SpikeDataNoLaserAll{v,j},numel(NOLASER_SPIKE_ALL{v,j}),0.001,[-Win Win],5);

    end
end


%Calculate Sparseness index for each cell


for v = 1:numcell
   for j = 1:nFreq
       mFR_ON(j) = mean(FR_LASER_ALL(j,v,241:290));
       mFR_OFF(j) = mean(FR_NOLASER_ALL(j,v,241:290));
   end

    LASERSidxON(v) = Sparseness(mFR_ON,nFreq);  
    LASERSidxOFF(v) = Sparseness(mFR_OFF,nFreq);  
end




figure; 
    subplot(1,2,1);
    scatter(LASERSidxOFF, LASERSidxON, 25,'filled')
    hold on; line([0 1], [0 1],'Color','k','LineStyle','--')
    hold off; axis square
    xlabel('Sparseness OFF')
    ylabel('Sparseness ON')
    set(gca,'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1])
    [p_LASERsparse, h_LASERsparse] = signrank(LASERSidxOFF, LASERSidxON);

    subplot(1,2,2);
    bar([nanmean(LASERSidxOFF) nanmean(LASERSidxON)],0.5,'EdgeColor','none');
    hold on; errorbar([nanmean(LASERSidxOFF) nanmean(LASERSidxON)],[nanstd(LASERSidxOFF)...
        ./sqrt(length(LASERSidxOFF)) nanstd(LASERSidxON)./sqrt(length(LASERSidxON))],'k','LineStyle','none','LineWidth',2)
    box off
    ylabel('Sparseness')
    set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out','YLim',[0 1]); axis square


set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 400]);
set(gcf,'Position',[0 0 800 400]);
suptitle([TITLE ': Affected cells'])
cd(FigOutput)
print(['LASERsparseness_' TITLE],'-dpdf','-r400')

cd(FileOutput)
save(TITLE,'h_LASERsparse','p_LASERsparse','LASERSidxON','LASERSidxOFF','-append')

disp('*******************************************************************');
disp(['8. Sparseness of ' num2str(size(h1,1)) ' affected cells: h = ' num2str(h_sparse) ' and p = ' num2str(p_sparse)]);
disp('*******************************************************************');


