%% Mouse and session numbers for different conditions (uncomment condition of interest before running).
clear;
% % SET 1: PVs, IC sessions ONLY (ChR2)
% TITLE = 'PV + ChR2 in IC (anesthetized)';
% MNum = [3100 3101 3102 3103 3104 3111 3112 3113 3121 3122 3124];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1,2]; Sesh{2} = [1,2]; Sesh{3} = [1,2]; Sesh{4} = [1,2]; Sesh{5} = [1,2];
% Sesh{6} = 1; Sesh{7} = [1,2]; Sesh{8} = 1; Sesh{9} = [2,3]; Sesh{10} = [2,3];
% Sesh{11} = 1;

% SET 2: PVs, AC sessions ONLY (ChR2)
% TITLE = 'PV + ChR2 in AC (anesthetized)';
% MNum = [3120 3121 3122 3124];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1]; Sesh{2} = [1]; Sesh{3} = [1]; Sesh{4} = [2];

% SET 3: SOMs, IC sessions ONLY (ChR2)
% TITLE = 'SOM + ChR2 in IC (anesthetized)';
% MNum = [3125 3126 3127 3135 3136 3137];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [2,3]; Sesh{2} = [1,2]; Sesh{3} = [2,3,4]; Sesh{4} = [2 3];
% Sesh{5} = [2 3]; Sesh{6} = [1];
 
% % SET 4: SOMs, AC sessions ONLY (ChR2)
% TITLE = 'SOM + ChR2 in AC (anesthetized)';
% MNum = [3127 3135 3136 3158 3159];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1; Sesh{4} = [1,2]; Sesh{5} = 1;
 
% % SET 5: CaMKIIs, IC sessions ONLY (ChR2)
% TITLE = 'CaMKII + ChR2 in IC (anesthetized)';
% MNum = [3128 3129 3130 3140 3141];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [2,3,4]; Sesh{2} = [2,3,4]; Sesh{3} = [2,3,4]; Sesh{4} = [2,3];
% Sesh{5} = [2,3];
 
% % SET 6: CaMKIIs, AC sessions ONLY (ChR2)
% TITLE = 'CaMKII + ChR2 in AC (anesthetized)';
% MNum = [3128 3129 3130 3140 3141];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1; Sesh{4} = 1; Sesh{5} = 1;

%%%%%%%%%%%%%
%   AWAKE
%%%%%%%%%%%%%

% SET 1: SOMs, IC only (ChR2)
% TITLE = 'SOM + ChR2 in IC (awake)';
% MNum = [3138 3145 3146 3176 3177];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1,2]; Sesh{2} = [1]; Sesh{3} = [1,2]; Sesh{4} = [1,2];
% Sesh{5} = [1,2];

% % SET 2: SOMs, AC only (ChR2)
% TITLE = 'SOM + ChR2 in AC (awake)';
% MNum = [3176 3177];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [3]; Sesh{2} = [3,4];

% % SET 3: CaMKIIs, IC only (ChR2)
% TITLE = 'CaMKII + ChR2 in IC (awake)';
% MNum = [3143 3150];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1]; Sesh{2} = [1,2];

% % SET 4: CaMKIIs, AC only (ChR2)
% TITLE = 'CaMKII + ChR2 in AC (awake)';
% MNum = [];
% Sesh = cell(1,length(MNum));


% % SET 5: PVs, IC only (ChR2)
% TITLE = 'PV + ChR2 in IC (awake)';
% MNum = [3149 3174 3175];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1,2]; Sesh{2} = [1 2]; Sesh{3} = [1 2];

% % SET 6: PVs, AC only (ChR2)
% TITLE = 'PV + ChR2 in AC (awake)';
% MNum = [3174 3175];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [3]; Sesh{2} = [3,4];

%% ArchT recordings
clear;
% % SET 1: PVs, IC sessions ONLY (ArchT)

% % SET 2: PVs, AC sessions ONLY (ArchT)

% % SET 3: SOMs IC sessions ONLY (ArchT)
% TITLE = 'SOM + ArchT in IC (anesthetized)';
% MNum = [3139 3144 3163];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [2,3]; Sesh{2} = [2,3]; Sesh{3} = 2;

% % SET 4: SOMs AC sessions ONLY (ArchT)
% TITLE = 'SOM + ArchT in AC (anesthetized)';
% MNum = [3139 3144 3160 3163];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [1]; Sesh{2} = [1]; Sesh{3} = [1,2]; Sesh{4} = 1;

% % SET 5: CaMKIIs, IC sessions ONLY (ArchT)
% TITLE = 'CaMKII + ArchT in IC (anesthetized)';
% MNum = [3154];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 4;

% % SET 6: CaMKIIs, AC sessions ONLY (ArchT)
% TITLE = 'CaMKII + ArchT in AC (anesthetized)';
% MNum = [3154];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [2,3];

%%%%%%%%%
% AWAKE
%%%%%%%%%

% % SET 1: PVs, IC only (ArchT)
% TITLE = 'PV + ArchT in IC (awake)';
% MNum = [3172 3182 3183];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 1; Sesh{2} = [1 2]; Sesh{3} = [1 2];

% % SET 2: PVs, AC only (ArchT)
% TITLE = 'PV + ArchT in AC (awake)';
% MNum = [3172 3182 3183];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 2; Sesh{2} = [3,4]; Sesh{3} = 3;

% % SET 3: SOMs, IC only (ArchT)
% TITLE = 'SOM + ArchT in IC (awake)';
% MNum = [3162 3164];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 1; Sesh{2} = [1,2];

% % SET 4: SOMs, AC only (ArchT)
TITLE = 'SOM + ArchT in AC (awake)';
MNum = [3161 3162];
Sesh = cell(1,length(MNum));
Sesh{1} = [1,2]; Sesh{2} = 2;

% % SET 5: CaMKIIs, IC only (ArchT)
% TITLE = 'CaMKII + ArchT in IC (awake)';
% MNum = [3153 3180 3181];
% Sesh = cell(1,length(MNum));
% Sesh{1} = 2; Sesh{2} = [1,2]; Sesh{3} = [1,2];

% % SET 6: CaMKIIs, AC only (ArchT)
% TITLE = 'CaMKII + ArchT in AC (awake)';
% MNum = [3180 3181];
% Sesh = cell(1,length(MNum));
% Sesh{1} = [3,4]; Sesh{2} = 3;

%% Calculate smoothed firing rate for each cell over highest amplitudes and top preferred frequencies

%Pre-allocate variables
allCELL = cell(1,length(MNum));
FR_LASER = cell(1,length(MNum)); 
FR_NOLASER = cell(1,length(MNum));
LOCSoff = cell(1,length(MNum));
LOCSon = cell(1,length(MNum));
CellQ = cell(1,length(MNum));
all_LASERSPIKE2 = cell(1,length(MNum)); 
all_NOLASERSPIKE2 = cell(1,length(MNum));


for u = 1:length(MNum)  
    %Pull out data files for mice and session numbers of interest
    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    h1 = ls('data\TC3-1*_laser.mat');
    idxUSE = find(ismember(str2num(h1(:,15)),Sesh{u}));
    h = h1(idxUSE,:);
    
    
    %Step 1: Use top 3 amplitudes, but must find the frequencies to use.   
    nfreq = 7; %Number of frequencies to use in smoothes response
    [LOCSon{u},LOCSoff{u}] = TC_Select_noGauss(MNum(u),h,nfreq,0);
    
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

        %Select those in bins around top frequencies and top 3 amplitudes
        LASER_SPIKE{v} = SpkTime_Laser(LOCSon{u}(v,:),6:8,:);
        NOLASER_SPIKE{v} = SpkTime_NoLaser(LOCSoff{u}(v,:),6:8,:);
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

cd('C:\Users\Jennifer\Documents\MATLAB\AuditoryCortex2017');
save(TITLE,'TITLE','MNum','Sesh','allCELL','FR_LASER','FR_NOLASER','CellQ');
%% Select good cells and calculate spontaneous firing rate and tone-evoked response magnitude for each cell

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

for u = 1:length(MNum)
    %Select good cells
    numcellA = size(allCELL{u},1);
    idxGOOD{u} = NaN(1,numcellA);    
    for j = 1:numcellA
        if CellQ{u}(j) <= 4 && ... %Use only good multi-units and single units
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

% figure; 
% subplot(1,2,1); scatter([spontOFF{:}], [spontON{:}],'*')
% hold on; line([0 max([[spontOFF{:}], [spontON{:}]])], [0 max([[spontOFF{:}], [spontON{:}]])]);
% hold off;
% xlabel('LASER OFF')
% ylabel('LASER ON')
% title('Spontaneous activity (Hz)')
% 
% subplot(1,2,2); scatter([magOFF{:}], [magON{:}],'*')
% hold on; line([0 max([[magOFF{:}], [magON{:}]])], [0 max([[magOFF{:}], [magON{:}]])]);
% hold off;
% xlabel('LASER OFF')
% ylabel('LASER ON')
% title('Tone-evoked response magnitude (Hz)')

NumberOfCells = length([IDX{:}])



%% Categorize each neuron based on spontaneous activity (to remove units with ONLY laser onset/offest responses)

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

% %Plot cells to check if categorization is working as expected.
% a = 3; b = 3; %Subplot dimensions (rows vs columns)
% nplot = a*b; %number of plots per figure
% nfig = ceil(nCell./nplot);
% allON = vertcat(FRon{:});
% allOFF = vertcat(FRoff{:});
% for xx = 1:nfig
%     for j = 1:nCell
%         if j < (xx*nplot + 1) && j > (xx - 1)*nplot
%             figure(xx); 
%             subplot(a,b,j - (xx-1)*nplot);
%             plot(allOFF(j,:),'k')
%             hold on; plot(allON(j,:),'m')
%             alpha(0.2)
%             hold off;
%             title([GOODCELLall(j,10:21) '; Cat ' num2str(CellCat(j))]);
%             xlabel('Time (ms)'); ylabel('FR (spikes/s)')
%             box off
% 
%         end
% 
%     end
%      set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition',[0 0 16 12], 'PaperSize', [16 12], 'color', 'white');
%      %print('-djpeg','-r400', [num2str(MNum(n)) '_allstim-' num2str(xx)]);
% end
% 
% %Display number of units per category
% Cat0 = length(find(CellCat == 0))
% Cat1 = length(find(CellCat == 1))
% Cat2 = length(find(CellCat == 2))
% Cat3 = length(find(CellCat == 3))

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

NumberOfCellsUse = length(spontOFFall2)

%% Use category 0 and category 1 cells to plot tone-evoked response magnitude and spontaneous activity
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

%figure('rend','painters','pos',[10 400 2000 400]); 
figure;
subplot(1,4,1); scatter(spontOFFall2, spontONall2,25,'filled')
LIM = max([spontOFFall2, spontONall2]);
LIM = 40;
hold on; line([0 LIM], [0 LIM],'Color','k','LineStyle','--');
hold off;
set(gca,'TickDir','out'); axis square
set(gca,'xlim',[0 LIM],'ylim',[0 LIM])
xlabel('LASER OFF')
ylabel('LASER ON')
title('Spontaneous activity (Hz)')

subplot(1,4,2); bar([mean(spontOFFall2) mean(spontONall2)],0.5,'EdgeColor','none');
hold on; errorbar([mean(spontOFFall2) mean(spontONall2)],[std(spontOFFall2)./sqrt(length(spontOFFall2)) std(spontONall2)./sqrt(length(spontONall2))],'k','LineStyle','none')
box off
ylabel('Firing rate (Hz)')
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
set(gca,'YLim',[0 12])
[p_spont, h_spont] = signrank(spontOFFall2, spontONall2)

subplot(1,4,3); scatter(magOFFall2, magONall2,25,'filled')
LIM = max([magOFFall2, magONall2]);
LIM = 100;
hold on; line([0 LIM], [0 LIM], 'Color','k','LineStyle','--');
hold off;
set(gca,'TickDir','out'); axis square
set(gca,'xlim',[0 LIM],'ylim',[-20 LIM])
xlabel('LASER OFF')
ylabel('LASER ON')
title('Tone-evoked response magnitude (Hz)')

subplot(1,4,4); bar([mean(magOFFall2) mean(magONall2)],0.5,'EdgeColor','none');
hold on; errorbar([mean(magOFFall2) mean(magONall2)],[std(magOFFall2)./sqrt(length(magOFFall2)) std(magONall2)./sqrt(length(magONall2))],'k','LineStyle','none')
box off
ylabel('Firing rate (Hz)')
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
set(gca,'YLim',[0 30])
[p_mag, h_mag] = signrank(magOFFall2, magONall2)

suptitle(TITLE);

cd C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\Conferences\AuditoryCortex2017\Figures
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 400]);
set(gcf,'Position',[0 0 1600 400]);
print(TITLE,'-dpdf','-r400')

%% Sparseness
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

nFreq = 50;
numcell = size(GOODCELLall,1);
FR_LASER_ALL = [];
FR_NOLASER_ALL = [];

h1 = GOODCELLall;
afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
Win = 0.24; %Window (in seconds) around start of tone to look for spikes

SpkTime_LaserAll = cell(size(h1,1),nFreq);
SpkTime_NoLaserAll = cell(size(h1,1),nFreq);

for v = 1:size(h1,1) %for filter 3, v = 47 peak is at index 2

    match = strfind(h1(v,:), '_');
    q = match(end);
    load(['D:\Spikes\M' h1(v,10:13) '\SpikeMat\R407F'  h1(v,5:q-1) '.mat']); 

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
[p_sparse, h_sparse] = signrank(SidxOFF, SidxON)

subplot(1,2,2);
bar([nanmean(SidxOFF) nanmean(SidxON)],0.5,'EdgeColor','none');
hold on; errorbar([nanmean(SidxOFF) nanmean(SidxON)],[nanstd(SidxOFF)./sqrt(length(SidxOFF)) nanstd(SidxON)./sqrt(length(SidxON))],'k','LineStyle','none')
box off
ylabel('Sparseness')
set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out','YLim',[0 1]); axis square

suptitle(TITLE)
cd C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\Conferences\AuditoryCortex2017\Figures
set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[800 400]);
set(gcf,'Position',[0 0 800 400]);
print(['sparseness_' TITLE],'-dpdf','-r400')

