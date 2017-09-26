set(0,'DefaultTextFontname', 'Arial')
%%
%close all;

% % SET 1: PVs, IC sessions ONLY (ChR2)
% MNum = [3100 3101 3102 3103 3104 3111 3112 3113 3121 3122 3124];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1,2]; Sesh{2} = [1,2]; Sesh{3} = [1,2]; Sesh{4} = [1,2]; Sesh{5} = [1,2];
% Sesh{6} = 1; Sesh{7} = [1,2]; Sesh{8} = 1; Sesh{9} = [2,3]; Sesh{10} = [2,3];
% Sesh{11} = 1;

% SET 2: PVs, AC sessions ONLY (ChR2)
% MNum = [3120 3121 3122 3124];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1]; Sesh{2} = [1]; Sesh{3} = [1]; Sesh{4} = [2];

% SET 3: SOMs, IC sessions ONLY (ChR2)
% MNum = [3125 3126 3127 3135 3136 3137];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [2,3]; Sesh{2} = [1,2]; Sesh{3} = [2,3,4]; Sesh{4} = [2 3];
% Sesh{5} = [2 3]; Sesh{6} = [1];
 
% % SET 4: SOMs, AC sessions ONLY (ChR2)
% MNum = [3127 3135 3136];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1;
 
% % SET 5: CaMKIIs, IC sessions ONLY (ChR2)
% MNum = [3128 3129 3130 3140 3141];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [2,3,4]; Sesh{2} = [2,3,4]; Sesh{3} = [2,3,4]; Sesh{4} = [2,3];
% Sesh{5} = [2,3];
 
% % SET 6: CaMKIIs, AC sessions ONLY (ChR2)
% MNum = [3128 3129 3130 3140 3141];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1; Sesh{4} = 1; Sesh{5} = 1;

%%%%%%%%%%%%%
%   AWAKE
%%%%%%%%%%%%%

% SET 1: SOMs, IC only (ChR2)
MNum = [3138 3145 3146];
Sesh = cell(1,length(MNum));
Sesh{1} = [1,2]; Sesh{2} = [1]; Sesh{3} = [1,2];

% % SET 2: CaMKIIs, IC only (ChR2)
% MNum = [3143 3150];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1]; Sesh{2} = [1,2];

% % SET 2: PVs, IC only (ChR2)
% MNum = [3149];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1,2];

%%
GOODCELL = cell(1,length(MNum));
FR_LASER = cell(1,length(MNum)); FR_NOLASER = cell(1,length(MNum));
LOCS = cell(1,length(MNum));
CellQ = cell(1,length(MNum));
all_LASERSPIKE2 = cell(1,length(MNum)); all_NOLASERSPIKE2 = cell(1,length(MNum));
for u = 1:length(MNum)   
    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    h1 = ls('data\TC3-1*_laser.mat');
    idxUSE = find(ismember(str2num(h1(:,15)),Sesh{u}));
    h = h1(idxUSE,:);
    %Step 1: Use top 3 amplitudes, but must find the frequencies to use.    
    [GaussParams,LOCS{u}] = TC_Select5(MNum(u),h,0);
    err = find( LOCS{u}.laser < 5 | LOCS{u}.nolaser < 5 | LOCS{u}.laser > 46 | LOCS{u}.nolaser > 46 | isnan(LOCS{u}.laser));        
    
    %Remove "bad" cells (that have peak too close to edge)
    badcell = err;
    if ~isempty(badcell)
        h(badcell,:) = []; 
        LOCS{u}.laser(badcell) = []; LOCS{u}.nolaser(badcell) = [];
    end
    
    %Step 2: Separate spikes by tone trial (separated by trial, frequency,
    %amplitude, and laser vs no laser)
    afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
    Win = 0.24; %Window (in seconds) around start of tone to look for spikes

    %Separate out spikes based on if during laser trials vs no laser trials
    LASER_SPIKE = cell(1,size(h,1));
    NOLASER_SPIKE = cell(1,size(h,1));
    for v = 1:size(h,1) 
        match = strfind(h(v,:), '0');
        q = match(end);
        load(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat\R407F'  h(v,5:q) '.mat']); 
        CellQ{u}(v) = CellInfo(6);

        [SpkTime_Laser, SpkTime_NoLaser, ~, ~] = SpikeTime(afo5,SpikeData,nRep,Win);

        %Select those in bins around best frequency and top 3 amplitudes
        LASER_SPIKE{v} = SpkTime_Laser(LOCS{u}.laser(v) - 3:LOCS{u}.laser(v) + 3,6:8,:);
        NOLASER_SPIKE{v} = SpkTime_NoLaser(LOCS{u}.nolaser(v) - 3:LOCS{u}.nolaser(v) + 3,6:8,:);
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
        
    GOODCELL{u} = h;
end


%% Calculate spontaneous firing rate and tone-evoked response magnitude for each cell
idxGOOD = cell(1,length(MNum));
spontON = cell(1,length(MNum)); spontOFF = cell(1,length(MNum)); magON = cell(1,length(MNum)); 
magOFF = cell(1,length(MNum)); aveON = cell(1,length(MNum)); aveOFF = cell(1,length(MNum));
for n = 1:length(MNum)
    numcell = size(GOODCELL{n},1);
    idxGOOD{n} = NaN(1,numcell);
    spontON{n} = NaN(1,numcell); spontOFF{n} = NaN(1,numcell); magON{n} = NaN(1,numcell);
    magOFF{n} = NaN(1,numcell); aveON{n} = NaN(1,numcell); aveOFF{n} = NaN(1,numcell);
    for j = 1:numcell
        if CellQ{n}(j) <= 4 && ... %Use only good multi-units and single units
            mean(FR_NOLASER{n}(j,240:290)) > (2*std(FR_NOLASER{n}(j,145:235)) + mean(FR_NOLASER{n}(j,145:235))); %Only use cells with ave tone-evoked FR higher than 2 std's above spontaneous
            idxGOOD{n}(j) = 1; %Identify usable cells
            spontON{n}(j) = mean(FR_LASER{n}(j,145:235)); spontOFF{n}(j) = mean(FR_NOLASER{n}(j,145:235));
            aveON{n}(j) = mean(FR_LASER{n}(j,240:315)); aveOFF{n}(j) = mean(FR_NOLASER{n}(j,240:315));
            magON{n}(j) = aveON{n}(j) - spontON{n}(j); magOFF{n}(j) = aveOFF{n}(j) - spontOFF{n}(j);
            

        end
    end
end

figure; 
scatter([spontOFF{:}], [spontON{:}],'*')
hold on; line([0 max([[spontOFF{:}], [spontON{:}]])], [0 max([[spontOFF{:}], [spontON{:}]])]);
hold off;
xlabel('LASER OFF')
ylabel('LASER ON')
title('Spontaneous activity (Hz)')

figure; 
scatter([magOFF{:}], [magON{:}],'*')
hold on; line([0 max([[magOFF{:}], [magON{:}]])], [0 max([[magOFF{:}], [magON{:}]])]);
hold off;
xlabel('LASER OFF')
ylabel('LASER ON')
title('Tone-evoked response magnitude (Hz)')

%% Choose units that are significantly affected by laser by a trial by trial comparison for each cell


%Identify neurons that show laser effect in each of 3 categories from
%previous section
FRon = cell(1,length(MNum));
FRoff = cell(1,length(MNum));
baseVAR = cell(1,length(MNum));
for n = 1:length(MNum)
    idx = find(~isnan(idxGOOD{n}));
    FRon{n} = FR_LASER{n}(idx,:);
    FRoff{n} = FR_NOLASER{n}(idx,:);
    baseVAR{n} = std(FRoff{n}(:,145:235),0,2); %baseline variability
end
Thresh = 1; %Threshold
baseVARall = vertcat(baseVAR{:})';

spontOFFall = [spontOFF{:}];
spontOFFall(isnan(spontOFFall)) = [];
spontONall = [spontON{:}];
spontONall(isnan(spontONall)) = [];
spontRatio1 = spontONall - (spontOFFall + Thresh*baseVARall);
spontRatio2 = spontONall - (spontOFFall - Thresh*baseVARall);
SigCellspontUP = find(spontRatio1 > 0); %Decrease activity with laser ON
SigCellspontDOWN = find(spontRatio2 < 0); %Increase activity with laser ON

magOFFall = [magOFF{:}];
magOFFall(isnan(magOFFall)) = [];
magONall = [magON{:}];
magONall(isnan(magONall)) = [];
magRatio1 = magONall - (magOFFall + Thresh*baseVARall);
magRatio2 = magONall - (magOFFall - Thresh*baseVARall);
SigCellmagUP = find(magRatio1 > 0);
SigCellmagDOWN = find(magRatio2 < 0);

%Re-plot ratios to see which have been identified as significant.
figure; 
subplot(3,2,1);
scatter(spontOFFall, spontONall,25,'filled')
%hold on; line([0 max([[spontOFF{:}], [spontON{:}]])], [0 max([[spontOFF{:}], [spontON{:}]])]);
hold on; line([0 31.98], [0 31.98]);
scatter(spontOFFall(SigCellspontDOWN), spontONall(SigCellspontDOWN),25,'m','filled')
scatter(spontOFFall(SigCellspontUP), spontONall(SigCellspontUP),25,'y','filled')
hold off; axis square
set(gca, 'xlim', [0 inf], 'ylim', [0 inf],'TickDir','out')
xlabel('LASER OFF'); ylabel('LASER ON')
title(['Spontaneous activity (Hz)'])
box off

subplot(3,2,2);
scatter(magOFFall, magONall,25,'filled')
%hold on; line([0 max([magOFFall, magONall])], [0 max([magOFFall, magONall])]);
hold on; line([0 63.56], [0 63.56]);
scatter(magOFFall(SigCellmagDOWN), magONall(SigCellmagDOWN),25,'m','filled')
scatter(magOFFall(SigCellmagUP), magONall(SigCellmagUP),25,'y','filled')
hold off; axis square
set(gca, 'xlim', [0 inf], 'ylim', [0 inf],'TickDir','out')
xlabel('LASER OFF')
title('Tone-evoked response magnitude (Hz)')
box off;


%Normalize FR:
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

%Plot average time-course for neurons identified as showing significant laser effect
subplot(3,2,3);
plot(mean(FR_NoLaser_M_Norm(SigCellspontDOWN,:),1),'k')
hold on; plot(mean(FR_Laser_M_Norm(SigCellspontDOWN,:),1),'m'); hold off;
ylabel('Normalized FR'); title(['N = ' num2str(length(SigCellspontDOWN))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:50:480,'XTickLabel',[-0.1:0.05:0.3])
subplot(3,2,5);
plot(mean(FR_NoLaser_M_Norm(SigCellspontUP,:),1),'k')
hold on; plot(mean(FR_Laser_M_Norm(SigCellspontUP,:),1),'y'); hold off;
ylabel('Normalized FR'); xlabel('Time (s)'); title(['N = ' num2str(length(SigCellspontUP))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:50:480,'XTickLabel',[-0.1:0.05:0.3])

subplot(3,2,4);
plot(mean(FR_NoLaser_M_Norm(SigCellmagDOWN,:),1),'k')
hold on; plot(mean(FR_Laser_M_Norm(SigCellmagDOWN,:),1),'m'); hold off;
title(['N = ' num2str(length(SigCellmagDOWN))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:50:480,'XTickLabel',[-0.1:0.05:0.3])
subplot(3,2,6);  
plot(mean(FR_NoLaser_M_Norm(SigCellmagUP,:),1),'k')
hold on; plot(mean(FR_Laser_M_Norm(SigCellmagUP,:),1),'y'); hold off;
xlabel('Time (s)'); title(['N = ' num2str(length(SigCellmagUP))]);
box off; set(gca,'TickDir','out','XLim',[0 500],'XTick',40:50:480,'XTickLabel',[-0.1:0.05:0.3])


%% Sparseness

nFreq = 50;
for u = 1:length(MNum)
    h1 = GOODCELL{u}(find(~isnan(idxGOOD{u})),:);
    afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
    Win = 0.24; %Window (in seconds) around start of tone to look for spikes

    SpkTime_LaserAll = cell(size(h1,1),nFreq);
    SpkTime_NoLaserAll = cell(size(h1,1),nFreq);

    for v = 1:size(h1,1) %for filter 3, v = 47 peak is at index 2

    match = strfind(h1(v,:), '0');
    q = match(end);
    load(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat\R407F'  h1(v,5:q) '.mat']); 

            [SpkTime_Laser, SpkTime_NoLaser, ~, ~] = SpikeTime(afo5,SpikeData,nRep,Win);
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
            FR_LASER_ALL{u}(j,v,:) = smoothFRx4(SpikeDataLaserAll{v,j},numel(LASER_SPIKE_ALL{v,j}),0.001,[-Win Win],5);

            %For no laser trials
            for w = 1:numel(SpkTime_NoLaserAll{v})
                NOLASER_SPIKE_ALL{v,j}{w} = [sort(SpkTime_NoLaserAll{v,j}{w}); w*ones(1,length(SpkTime_NoLaserAll{v,j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                SpikeDataNoLaserAll{v,j} = [SpikeDataNoLaserAll{v,j} NOLASER_SPIKE_ALL{v,j}{w}];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
            end
            FR_NOLASER_ALL{u}(j,v,:) = smoothFRx4(SpikeDataNoLaserAll{v,j},numel(NOLASER_SPIKE_ALL{v,j}),0.001,[-Win Win],5);

        end
    end
   

end

FR_LaserAll_M = [FR_LASER_ALL{:}];
FR_NoLaserAll_M = [FR_NOLASER_ALL{:}];

%Calculate Sparseness index for each cell
numcell = nansum([idxGOOD{:}]);
SidxON = nan(1,numcell);
SidxOFF = nan(1,numcell);
for v = 1:numcell
   for j = 1:nFreq
       mFR_ON(j) = mean(FR_LaserAll_M(j,v,241:290));
       mFR_OFF(j) = mean(FR_NoLaserAll_M(j,v,241:290));
   end
       
    SidxON(v) = Sparseness(mFR_ON,nFreq);  
    SidxOFF(v) = Sparseness(mFR_OFF,nFreq);  
end

figure; subplot(1,2,1);
scatter(SidxOFF, SidxON, '.')
hold on; line([0 1], [0 1])
hold off; axis square
xlabel('Sparseness OFF')
ylabel('Sparseness ON')
title('Sparseness')
set(gca,'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1])
[p_sparse, h_sparse] = signrank(SidxOFF, SidxON)

subplot(1,2,2);
bar([0.5 0.6], [nanmean(SidxOFF) nanmean(SidxON)],0.2);
hold on; errorbar([0.5,0.6],[nanmean(SidxOFF) nanmean(SidxON)],[nanstd(SidxOFF)./sqrt(length(SidxOFF)) nanstd(SidxON)./sqrt(length(SidxON))],'k.')
box off
ylabel('Sparseness')
set(gca,'XTick',[0.5 0.6], 'XTickLabels',{'OFF';'ON'},'TickDir','out','XLim',[0.47 0.63],'YLim',[0 1])

%% Plot tuning indices for FM sweeps (already calculated during initial analysis using FMsweepTC and FMsweepTC2)

% % SET 1: PVs, IC sessions ONLY (ChR2)
% MNum = [3102 3103 3104 3111 3112 3113 3121 3122 3124];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1,2]; Sesh{2} = [1,2]; Sesh{3} = [1,2];
% Sesh{4} = 1; Sesh{5} = [1,2]; Sesh{6} = 1; Sesh{7} = [2,3]; Sesh{8} = [2,3];
% Sesh{9} = 1;

TUNEIDX = cell(1,length(MNum));
DIRIDX = cell(1,length(MNum));
for u = 1:length(MNum)
    cd(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat']);
    FileList_temp = ls('J001*');
    idxUSE = find(ismember(str2num(FileList_temp(:,16)),Sesh{u}));
    FileList = FileList_temp(idxUSE,:);
    
    TUNEIDX{u} = NaN(size(FileList,1),2);
    DIRIDX{u} = NaN(size(FileList,1),2);
    for v = 1:size(FileList,1)
        load( FileList(v,:) );
        if CellInfo(6) <= 4 %&& (SWEEPS.tuneIDX(1) > mean(SWEEPS.bootTUNE(:,1)) + 0.5*std(SWEEPS.bootTUNE(:,1)) ...
                %|| SWEEPS.tuneIDX(1) > mean(SWEEPS.bootTUNE(:,1)) - 0.5*std(SWEEPS.bootTUNE(:,1))) 
           TUNEIDX{u}(v,:) = SWEEPS.tuneIDX;           
        end
        if CellInfo(6) <= 4 %&& (SWEEPS.dirIDX(1) > mean(SWEEPS.bootDIR(:,1)) + 0.5*std(SWEEPS.bootDIR(:,1)) ...
               % || SWEEPS.dirIDX(1) > mean(SWEEPS.bootDIR(:,1)) - 0.5*std(SWEEPS.bootDIR(:,1))) 
           DIRIDX{u}(v,:) = SWEEPS.dirIDX;          
        end
        
    end
    
end

TUNEIDXall = vertcat(TUNEIDX{:});
DIRIDXall = vertcat(DIRIDX{:});

figure;
subplot(1,2,1);
scatter(TUNEIDXall(:,1),TUNEIDXall(:,2),'.')
hold on; line([0 1],[0 1]); hold off;
xlabel('LASER OFF'); ylabel('LASER ON'); title('SPARSENESS')
axis square; set(gca,'TickDir','out', 'YTick',[0 0.5 1])
[p_tune, h_tune] = signrank(TUNEIDXall(:,1),TUNEIDXall(:,2))

subplot(1,2,2);
scatter(DIRIDXall(:,1),DIRIDXall(:,2),'.')
hold on; line([-1 1],[-1 1]); hold off;
xlabel('LASER OFF'); ylabel('LASER ON'); title('DIRECTION INDEX')
axis square; set(gca,'TickDir','out')
[p_dir, h_dir] = signrank(DIRIDXall(:,1) - DIRIDXall(:,2))



