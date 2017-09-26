


%% Select good cells (based on one-sided t-test comparing sequential laser/nolaser trials from 400-450ms [spontaneous firing])
clear all

MNum = [3080];
Filters = [1 2]; %Session 1 = 290; Session 2 = 110; Session 3 = 60;

SIG = cell(length(MNum),1);
for u = 1:length(MNum)

    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    filnam = ls('data\TC3-1*_laser.mat');
    filecut = filnam(:,5:end-12);
    h1 = unique(filecut,'rows');

    %Now calculate spike times for laser vs no laser trials
    for v = 1:size(h1,1) 
            spkcountNL = [];
            spkcountL = [];
            spkcount = [];
        for i = 1:length(Filters)
            %filt = str2num(h1(v,11));
            match = strfind(h1(v,:), '0');
            q = match(end);
            load(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat\R407F'  h1(v,1:q) '_' num2str(Filters(i)) '.mat']); 
            CellQ(u,v) = CellInfo(end); %Save cell quality

            %End times for each trial
            Time = [0.5:0.5:nRep*400]; 
            T_Laser = Time(1:2:nRep*400);
            T_NoLaser = Time(2:2:nRep*400);
            

            %For laser and no laser trials, count number of spikes from
            %first 50ms of when laser would be on

            for ii = 1:length(T_Laser)
               spkcountNL(ii,i) = numel(find(SpikeData(2,:) < (T_NoLaser(ii)-0.05) & SpikeData(2,:) > (T_NoLaser(ii) - 0.1))); 
               spkcountL(ii,i) = numel(find(SpikeData(2,:) < (T_Laser(ii)-0.05) & SpikeData(2,:) > (T_Laser(ii) - 0.1)));                
            end
            

        end
        %Include all 3 LED powers
        spkcount(:,1) = [spkcountNL(:,1); spkcountNL(:,2)]; 
        spkcount(:,2) = [spkcountL(:,1); spkcountL(:,2)];
        
        %spkcount(:,1) = [spkcountNL(:,1)]; 
        %spkcount(:,2) = [spkcountL(:,1)];

        [h,p] = ttest(spkcount(:,1),spkcount(:,2),'Tail','right'); %h = 1 indicates no laser is LARGER than laser
        SIG{u}(v,:) = [h,p];
    end

    %Now, select cells based on effect of laser on spontaneous FR
    G = zeros(1,length(h1));
    for v = 1:length(h1)
        if SIG{u}(v,1) == 1 && CellQ(u,v) < 5 
            G(v) = 1;
        end
    end
    numGoodCell = sum(G)

    idx_cell{u} = find(G == 1);
    gdcell{u} = h1(idx_cell{u},:);
        
end

%FOR M3082, FOR FILTER 3 ALONE THERE ARE 6 SIGNIFICANT CELLS ONLY, of those
%only 3 are usable ( 3 give errors due to TC peak location )
%FOR M3082, FOR FILTER 1,2 COMBINED THERE ARE 32 SIGNIFICANT CELLS

%% Compute the difference between tone-evoked and baseline firing rate (separate for each laser power)
% MNum = [3015 3016 3042 3043 3044 3050 3053]; 
 %Filters = [1];

for u = 1:length(MNum)%1:length(MNum)     
    %h = GOODCELLS{u};
    h = gdcell{u};
    if ~isempty(h)
%Step 1: Use top 3 amplitudes, but must find the frequencies to use.    
    for i = 1:length(Filters)
        [GaussParams{i},LOCS{i}] = TC_Select3(MNum(u),h,Filters(i),0);
        err{i} = find( LOCS{i}.laser < 5 | LOCS{i}.nolaser < 5 | LOCS{i}.laser > 46 | LOCS{i}.nolaser > 46 | isnan(LOCS{i}.laser));        
    end

    % If bad for one filter, remove for all of them.
    badcell = unique([err{:}]);
    if ~isempty(badcell)
        h(badcell,:) = []; 
        for i = 1:length(Filters)
            LOCS{i}.laser(badcell) = []; LOCS{i}.nolaser(badcell) = [];
            GaussParams{i}.rsqON(badcell) = []; GaussParams{i}.rsqOFF(badcell) = [];
            GaussParams{i}.BestFreqOFF(badcell) = []; GaussParams{i}.BestFreqON(badcell) = [];
            GaussParams{i}.WidthON(badcell) = []; GaussParams{i}.WidthOFF(badcell) = [];
        end
    end
 
    
    %Step 2: Separate spikes by tone trial (separated by trial, frequency,
    %amplitude, and laser vs no laser)
    afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
    Win = 0.24; %Window (in seconds) around start of tone to look for spikes


    %Separate out spikes based on if during laser trials vs no laser trials
    for i = 1:length(Filters)
        for v = 1:size(h,1) 
            match = strfind(h(v,5:end), '0');
            q = match(end);
            load(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat\R407F'  h(v,1:5+q-1) '_' num2str(Filters(i)) '.mat']); 
            %CellQ(u,v) = CellInfo(end); %Save cell quality

            [SpkTime_Laser{i}, SpkTime_NoLaser{i}, ~, ~] = SpikeTime(afo5,SpikeData,nRep,Win);

            %Select those in bins around best frequency and top 3 amplitudes
            LASER_SPIKE{v,i} = SpkTime_Laser{i}(LOCS{i}.laser(v) - 3:LOCS{i}.laser(v) + 3,6:8,:);
            NOLASER_SPIKE{v,i} = SpkTime_NoLaser{i}(LOCS{i}.nolaser(v) - 3:LOCS{i}.nolaser(v) + 3,6:8,:);
        end
    end

    %Step 3: Calculate and smooth firing rate for laser and no laser conditions
    SpikeDataLaser = cell(1,length(h));        
    SpikeDataNoLaser = cell(1,length(h));

    for i = 1:length(Filters)
        clear LASER_SPIKE2
        for v = 1:size(h,1)
            for w = 1:numel(LASER_SPIKE{v,i})
                LASER_SPIKE2{v}{w} = [sort(LASER_SPIKE{v,i}{w}); w*ones(1,length(LASER_SPIKE{v,i}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                SpikeDataLaser{v} = [SpikeDataLaser{v} LASER_SPIKE2{v}{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
            end
            FR_LASER{u,i}(v,:) = smoothFRx2(SpikeDataLaser{v},numel(LASER_SPIKE{v,i}),0.001,[-Win Win],[50,10]);
            
            for w = 1:numel(NOLASER_SPIKE{v,i})
                LASER_SPIKE2{v}{w} = [sort(NOLASER_SPIKE{v,i}{w}); w*ones(1,length(NOLASER_SPIKE{v,i}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                SpikeDataNoLaser{v} = [SpikeDataNoLaser{v} LASER_SPIKE2{v}{w}]; %Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
            end
            FR_NOLASER{u,i}(v,:) = smoothFRx2(SpikeDataNoLaser{v},numel(NOLASER_SPIKE{v,i}),0.001,[-Win, Win],[50,10]);
        end
    end
    
    %Average Laser and No Laser trials across conditions filter = 0 and filter = 7    
    for w = 1:size(h,1)
         FR_Laser_Mn{u}(w,:) = mean([FR_LASER{u,1}(w,:); FR_LASER{u,2}(w,:)]);
         FR_NoLaser_Mn{u}(w,:) = mean([FR_NOLASER{u,1}(w,:); FR_NOLASER{u,2}(w,:)]);
        
        %FR_Laser_Mn{u}(w,:) = FR_LASER{u,1}(w,:);
        %FR_NoLaser_Mn{u}(w,:) =FR_NOLASER{u,1}(w,:);
    end

    %Selection criterion #2: Mean firing rate during tone presentation must
    %be (thresh2 value) standard deviations above baseline firing rate
    G = zeros(1,size(h,1));
    thresh2 = 1.8;
    for v = 1:size(h,1)
        MN1 = mean(FR_NoLaser_Mn{u}(v,160:220));
        SD = std(FR_NoLaser_Mn{u}(v,160:220));
        MN2 = mean(FR_NoLaser_Mn{u}(v,241:290));
        if MN2 > MN1 + thresh2*SD %
            G(v) = 1;
        end
    end


        
    number_of_good_cells = sum(G) %Display number of good cells found (not necessary, but nice to see)
    ind_cell{u} = find(G == 1); %Adjust indices so we know which ones out of this data to plot
    GOODCELLS{u} = h(ind_cell{u},:);
    
    

    %Pull out best frequency calculated by TC_Select3 :
    for i = 1:length(Filters)
        for v = 1:size(h,1)
            BestFreq.laser{u}(i,v) = GaussParams{i}.BestFreqON(v);
            BestFreq.nolaser{u}(i,v) = GaussParams{i}.BestFreqOFF(v);
        end
    end
    
    %Pull out tuning curve width calculated by TC_Select3
    for i = 1:length(Filters)
        for v = 1:size(h,1)
            Width.laser{u}(i,v) = GaussParams{i}.WidthON(v);
            Width.nolaser{u}(i,v) = GaussParams{i}.WidthOFF(v);
        end
    end    
    
    else GOODCELLS{u} = []; ind_cell{u} = [];
    end
end
%% Make plots
ii = 2;
% Filters = [0 7 50];
 for n = 1%1%1%2:7
     cd(['D:\Spikes\M' num2str(MNum(n)) '\'])
     close all
     numcell = length(ind_cell{n}); %Number of good cells over both animals
     
%      for i = 1:length(Filters)
%          FR_Laser_M{i} = [FR_LASER{n,i}(ind_cell{n},:)];
%          FR_NoLaser_M{i} = [FR_NOLASER{n,i}(ind_cell{n},:)];
%      end
% 
%      %Take mean of laser and no laser across all 3 conditions
%      FR_Laser_Mn = [];
%      FR_NoLaser_Mn = [];
%      for w = 1:length(ind_cell{n})
%         FR_Laser_Mn(w,:) = mean([FR_Laser_M{1}(w,:); FR_Laser_M{2}(w,:)]);
%         FR_NoLaser_Mn(w,:) = mean([FR_NoLaser_M{1}(w,:); FR_NoLaser_M{2}(w,:); FR_NoLaser_M{3}(w,:)]);
%      end
    
    a = 3; b = 3; %Subplot dimensions (rows vs columns)
    nplot = a*b; %number of plots per figure
    nfig = ceil(length(ind_cell{n})./nplot); %Max number of figures needed to plot all cells from mouse n

    for xx = 1:nfig
        for j = 1:numcell
            if j < (xx*nplot + 1) && j > (xx - 1)*nplot
                figure(xx); 
                subplot(a,b,j - (xx-1)*nplot);
                plot(FR_NOLASER{n,ii}(j,:),'k')
                hold on; plot(FR_LASER{n,ii}(j,:),'m')
                MX = max([FR_NOLASER{n,ii}(j,:) FR_LASER{n,ii}(j,:)]); %Find max FR out of both laser and no laser
                %area([140 390],[MX MX],'LineStyle','none'); %Location of laser on
                %area([240 290],[MX MX],'LineStyle','none','FaceColor','r'); %Location of tone on
                alpha(0.2)
                hold off;
                title(['Cell ' num2str(j)]);
                xlabel('Time (ms)'); ylabel('FR (spikes/s)')
                box off
                
            end
            
        end
         set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition',[0 0 16 12], 'PaperSize', [16 12], 'color', 'white');
         %print('-djpeg','-r400', [num2str(MNum(n)) '_allstim-' num2str(xx)]);
     end
end
% %Plot average firing rate for each mouse across all cells 

 numcell = sum(cellfun('length',ind_cell)) %Calculate number of good cells across animals
% %Concatenate across animals 
for u = 1:length(MNum)
    FR_Laser_gd{u} = FR_Laser_Mn{u}(ind_cell{u},:)';
    FR_NoLaser_gd{u} = FR_NoLaser_Mn{u}(ind_cell{u},:)';
end
FR_Laser_M = [FR_Laser_gd{:}]';
FR_NoLaser_M = [FR_NoLaser_gd{:}]';

%Normalize firing rate
BaseFR = [];
ToneFR = [];
FR_Laser_M_Norm = [];
FR_NoLaser_M_Norm = [];
for j = 1:numcell
    BaseFR(j) = mean(FR_NoLaser_M(j, 181:230));
    ToneFR(j) = max(FR_NoLaser_M(j,241:290));
    
    FR_Laser_M_Norm(j,:) = (FR_Laser_M(j,:) - BaseFR(j))/(ToneFR(j) - BaseFR(j));
    FR_NoLaser_M_Norm(j,:) = (FR_NoLaser_M(j,:) - BaseFR(j))/(ToneFR(j) - BaseFR(j));
end


figure;
plot(mean(FR_NoLaser_M_Norm,1),'k')
hold on; plot(mean(FR_Laser_M_Norm,1),'m')
line([140 140], [-0.1 1],'LineStyle','--'); 
line([390 390], [-0.1 1],'LineStyle','--'); 
line([240 240],[-0.1 1],'LineStyle','--');
line([290 290],[-0.1 1],'LineStyle','--');
alpha(0.2)
hold off;
xlabel('Time'); ylabel('Normalized FR')
box off

%% Plot scatter (and bar) of Tone response
%Tone response (difference between baseline and peak tone FR)

ToneRespL = mean(FR_Laser_M_Norm(:,241:290), 2) - mean(FR_Laser_M_Norm(:,181:230), 2);
ToneRespNL = mean(FR_NoLaser_M_Norm(:,241:290), 2) - mean(FR_NoLaser_M_Norm(:,181:230), 2);

figure;
scatter(ToneRespNL,ToneRespL,'.')
hold on; line([0 1],[0 1])
axis square
hold off
xlabel('Tone response OFF')
ylabel('Tone response ON')
[p_tone, h_tone] = signrank(ToneRespNL, ToneRespL)


figure; bar([1 2],[mean(ToneRespNL) mean(ToneRespL)],0.6);
hold on; errorbar([1 2],[mean(ToneRespNL) mean(ToneRespL)],[std(ToneRespNL)./sqrt(length(ToneRespNL)) std(ToneRespL)./sqrt(length(ToneRespL))],'k.')
box off
ylabel('Tone Response')
set(gca,'XTickLabels',{'OFF';'ON'})


BaseRespL = mean(FR_Laser_M(:,181:230),2);
BaseRespNL = mean(FR_NoLaser_M(:,181:230),2);

figure;
scatter(BaseRespNL,BaseRespL,'.')
hold on; line([0 16],[0 16])
axis square
hold off
xlabel('Baseline FR OFF')
ylabel('Baseline FR ON')
[p_base,h_base] = signrank(BaseRespNL, BaseRespL)

figure; bar([1 2],[mean(BaseRespNL) mean(BaseRespL)],0.6);
hold on; errorbar([1 2],[mean(BaseRespNL) mean(BaseRespL)],[std(BaseRespNL)./sqrt(length(BaseRespNL)) std(BaseRespL)./sqrt(length(BaseRespL))],'k.')
box off
ylabel('Baseline activity (Hz)')
set(gca,'XTickLabels',{'OFF';'ON'})


%% Plot best frequency laser ON vs best frequency laser OFF (pool across laser conditions used in analysis)

for u = 1:length(MNum)
    BFuse.laser{u} = BestFreq.laser{u}(:,ind_cell{u});
    BFuse.nolaser{u} = BestFreq.nolaser{u}(:,ind_cell{u});
end

BFall.laser = [BFuse.laser{:}];
BFall.nolaser = [BFuse.nolaser{:}];

BFall.laser = reshape(BFall.laser',1,size(BFall.laser,1)*size(BFall.laser,2));
BFall.nolaser = reshape(BFall.nolaser',1,size(BFall.nolaser,1)*size(BFall.nolaser,2));


figure; scatter(BFall.nolaser,BFall.laser,'.')
hold on; line([0 max(BFall.nolaser, BFall.laser)],[0 max(BFall.nolaser, BFall.laser)])
axis square
xlabel('BF laser OFF (Hz)')
ylabel('BF laser ON (Hz)')
title('Best Frequency')
[p_BF, h_BF] = signrank(BFall.nolaser,BFall.laser)


%% Plot tuning curve width laser ON vs best frequency laser OFF (pool across laser conditions used in analysis)

for u = 1:length(MNum)
    WIDTHuse.laser{u} = Width.laser{u}(:,ind_cell{u});
    WIDTHuse.nolaser{u} = Width.nolaser{u}(:,ind_cell{u});
end

WIDTHall.laser = [WIDTHuse.laser{:}];
WIDTHall.nolaser = [WIDTHuse.nolaser{:}];

WIDTHall.laser = reshape(WIDTHall.laser',1,size(WIDTHall.laser,1)*size(WIDTHall.laser,2));
WIDTHall.nolaser = reshape(WIDTHall.nolaser',1,size(WIDTHall.nolaser,1)*size(WIDTHall.nolaser,2));

WIDTHall.laser(WIDTHall.laser < 0) = NaN;
WIDTHall.nolaser(WIDTHall.nolaser < 0) = NaN;


figure; scatter(WIDTHall.nolaser,WIDTHall.laser,'.')
hold on; line([0 max(WIDTHall.nolaser, WIDTHall.laser)],[0 max(WIDTHall.nolaser, WIDTHall.laser)])
axis square
xlabel('Tuning curve width OFF (Hz)')
ylabel('Tuning curve width ON (Hz)')
title('Tuning curve bandwidth')
[p_width, h_width] = ttest(WIDTHall.nolaser,WIDTHall.laser)

%% Sparseness
%load('C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\SfN2015\ttestGoodCells2')
nFreq = 50;
for u = 1:length(MNum)
    h1 = GOODCELLS{u};
    afo5 = load('D:\Code\TuningCurve\R407_pars'); %Load stimulus parameters
    Win = 0.24; %Window (in seconds) around start of tone to look for spikes

    SpkTime_LaserAll = cell(size(h1,1),nFreq);
    SpkTime_NoLaserAll = cell(size(h1,1),nFreq);
    for i = 1:length(Filters)
        for v = 1:size(h1,1) %for filter 3, v = 47 peak is at index 2

                match = strfind(h(v,5:end), '0');
                q = match(end);
                load(['D:\Spikes\M' num2str(MNum(u)) '\SpikeMat\R407F'  h(v,1:5+q-1) '_' num2str(Filters(i)) '.mat']); 

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
                FR_LASER_ALL{u,i}(j,v,:) = smoothFRx2(SpikeDataLaserAll{v,j},numel(LASER_SPIKE_ALL{v,j}),0.001,[-Win Win],[50,10]);

                %For no laser trials
                for w = 1:numel(SpkTime_NoLaserAll{v})
                    NOLASER_SPIKE_ALL{v,j}{w} = [sort(SpkTime_NoLaserAll{v,j}{w}); w*ones(1,length(SpkTime_NoLaserAll{v,j}{w}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                    SpikeDataNoLaserAll{v,j} = [SpikeDataNoLaserAll{v,j} NOLASER_SPIKE_ALL{v,j}{w}];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
                end
                FR_NOLASER_ALL{u,i}(j,v,:) = smoothFRx2(SpikeDataNoLaserAll{v,j},numel(NOLASER_SPIKE_ALL{v,j}),0.001,[-Win Win],[50,10]);

            end
        end
   
    end  
    for j = 1:nFreq
        for v = 1:size(h1,1)
            %FR_LaserAll_Mn{u}(j,v,:) = mean([FR_LASER_ALL{u,1}(j,v,:); FR_LASER_ALL{u,2}(j,v,:)]);
            %FR_NoLaserAll_Mn{u}(j,v,:) = mean([FR_NOLASER_ALL{u,1}(j,v,:); FR_NOLASER_ALL{u,2}(j,v,:)]);
            FR_LaserAll_Mn{u}(j,v,:) = FR_LASER_ALL{u,1}(j,v,:);
            FR_NoLaserAll_Mn{u}(j,v,:) = FR_NOLASER_ALL{u,1}(j,v,:);
        end
    end   
end

FR_LaserAll_M = [FR_LaserAll_Mn{:}];
FR_NoLaserAll_M = [FR_NoLaserAll_Mn{:}];

%Calculate Sparseness index for each cell
numcell = sum(cellfun('size',GOODCELLS,1));
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

figure; scatter(SidxOFF, SidxON, '.')
hold on; line([0 1], [0 1])
hold off; axis square
xlabel('Sparseness OFF')
ylabel('Sparseness ON')
title('Sparseness')
[h_sparse, p_sparse] = signrank(SidxOFF, SidxON)

figure; bar([1 2],[nanmean(SidxOFF) nanmean(SidxON)],0.6);
hold on; errorbar([1 2],[nanmean(SidxOFF) nanmean(SidxON)],[nanstd(SidxOFF)./sqrt(length(SidxOFF)) nanstd(SidxON)./sqrt(length(SidxON))],'k.')
box off
ylabel('Sparseness')
set(gca,'XTickLabels',{'OFF';'ON'})



%% AVERAGE TUNING CURVE: (QUESTION NEEDS TO BE ADDRESSED ABOUT NORMALIZATION)
%load('C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\SfN2015\ttestGoodCells2') %Load file names for good cells selected for analysis
FR1_laser = cell(length(MNum),length(Filters));
FR1_nolaser = cell(length(MNum),length(Filters));
for u = 1:length(MNum)
    h = GOODCELLS{u};
    if ~isempty(h)
        laser = cell(1,length(Filters));
        nolaser = cell(1,length(Filters));
        f1 = cell(1,length(Filters));
        for v = 1:size(h,1)
            match = strfind(h(v,5:end), '0');
            q = match(end);
            for i = 1:length(Filters)
                %Calculate smoothed tuning curves
                load(['D:\Spikes\M' num2str(MNum(u)) '\TCs\data\TC3-'  h(v,1:5+q-1) '_' num2str(Filters(i)) '_laser.mat']); 

                for w = 1:8
                    T1(:, w) = SmoothGaus(TC1.TCmat{1}(:, w), 3);
                end

                for w = 1:8
                    T2(:, w) = SmoothGaus(TC2.TCmat{1}(:, w), 3);
                end

                %Take mean across top 3 amplitudes
                laser{i}(:,v) = mean(T1(:,6:8),2);
                nolaser{i}(:,v) = mean(T2(:,6:8),2);
                
                %Step 1: Convert frequencies to octaves relative to center frequency
                fB_laser = BFuse.laser{u}(i,v); fB_nolaser = BFuse.nolaser{u}(i,v);
                f1_laser{i}(:,v) = log2(TC1.F./fB_laser); %Transform frequencies to be octaves relative to best frequency
                f1_nolaser{i}(:,v) = log2(TC1.F./fB_nolaser);
                
                %Step 2: Interpolate. [Interp 1 adds only 12 new points]
                f2 = [-3:0.1:3];
                FR1_laser{u,i}(:,v) = interp1(f1_laser{i}(:,v),laser{i}(:,v),f2);
                FR1_nolaser{u,i}(:,v) = interp1(f1_nolaser{i}(:,v),nolaser{i}(:,v),f2);
                
                %Normalize interpolated tuning curves: (QUESTION: Normalize
                %both by same number???)
                %normFR1_laser{u,i}(:,v) = (FR1_laser{u,i}(:,v) - min(FR1_laser{u,i}(:,v)))./(max(FR1_laser{u,i}(:,v))-min(FR1_laser{u,i}(:,v)));
                %normFR1_nolaser{u,i}(:,v) = (FR1_nolaser{u,i}(:,v)-min(FR1_nolaser{u,i}(:,v)))./(max(FR1_nolaser{u,i}(:,v))-min(FR1_nolaser{u,i}(:,v)));
                
            end
        end        

    end
end

%a = reshape(normFR1_laser,1,numel(normFR1_laser));
a = reshape(FR1_laser,1,numel(FR1_laser));
catFR1_laser = [a{:}];

%a = reshape(normFR1_nolaser,1,numel(normFR1_nolaser));
a = reshape(FR1_nolaser,1,numel(FR1_nolaser));
catFR1_nolaser = [a{:}];

meanFR1_laser = nanmean( catFR1_laser ,2 );
meanFR1_nolaser = nanmean(catFR1_nolaser,2);

meanFR1_laser_norm = (meanFR1_laser - min(meanFR1_laser))./(max(meanFR1_laser) - min(meanFR1_laser));
meanFR1_nolaser_norm = (meanFR1_nolaser - min(meanFR1_nolaser))./(max(meanFR1_nolaser) - min(meanFR1_nolaser));

figure; plot(f2,meanFR1_nolaser_norm,'k')
hold on; plot(f2,meanFR1_laser_norm,'r');
hold off; box off;
ylabel('Normalized FR')
xlabel('Octaves from best frequency')

%% AVERAGE TUNING CURVE (using sdev's from bf):
load('C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\SfN2015\ttestGoodCells2') %Load file names for good cells selected for analysis
FR1_laser = cell(length(MNum),length(Filters)-1);
FR1_nolaser = cell(length(MNum),length(Filters)-1);
for u = 1:length(MNum)
    h = GOODCELLS{u};
    if ~isempty(h)
        laser = cell(1,length(Filters) - 1);
        nolaser = cell(1,length(Filters) - 1);
        f1 = cell(1,length(Filters) - 1);
        for v = 1:size(h,1)
            match = strfind(h(v,5:end), '0');
            q = match(end);
            for i = 1:length(Filters) - 1
                %Calculate smoothed tuning curves
                load(['D:\Spikes\M' num2str(MNum(u)) '\TCs\data\' h(v,1:5+q-1) '-' num2str(Filters(i)) '_laser.mat']);

                for w = 1:8
                    T1(:, w) = SmoothGaus(TC1.TCmat{1}(:, w), 3);
                end

                for w = 1:8
                    T2(:, w) = SmoothGaus(TC2.TCmat{1}(:, w), 3);
                end

                %Take mean across top 3 amplitudes
                laser{i}(:,v) = mean(T1(:,6:8),2);
                nolaser{i}(:,v) = mean(T2(:,6:8),2);
                
                x = TC1.F';
                norm_func = @(b,x) b(1) .* exp(-((x - b(2))/b(3)).^2) + b(4);
                initguess = [5 15000 5000 0]; %Inital guess of paramters [amplitude, center freq, sdev, baseline shift]
                coeffval_LaserOFF = nlinfit(x,nolaser{i}(:,v),norm_func,initguess); 
                coeffval_LaserON = nlinfit(x,laser{i}(:,v),norm_func,initguess);

                %Start pulling out paramters:
                bf_ON = coeffval_LaserON(2); %Best frequency
                f2_ON = coeffval_LaserON(2)+coeffval_LaserON(3)*sqrt(log(2)); %higher frequensy at 1/2 max amplitude

                bf_OFF = coeffval_LaserOFF(2);
                f2_OFF = coeffval_LaserOFF(2)+coeffval_LaserOFF(3)*sqrt(log(2)); %higher frequency at 1/2 max amplitude
                
                %Step 1: Convert frequencies to octaves relative to center frequency
                f1_laser{i}(:,v) = log2(TC1.F./bf_ON); %Transform frequencies to be octaves relative to best frequency
                f1_nolaser{i}(:,v) = log2(TC1.F./bf_OFF);
                
                delf_laser = log2(f2_ON/bf_ON); delf_nolaser = log2(f2_OFF/bf_OFF);
                f2_laser{i}(:,v) = f1_laser{i}(:,v)./delf_laser;
                f2_nolaser{i}(:,v) = f1_nolaser{i}(:,v)./delf_nolaser;
                
                %Step 2: Interpolate. [Interp 1 adds only 12 new points]
                f3 = [-5:0.1:5];
                FR1_laser{u,i}(:,v) = interp1(f2_laser{i}(:,v),laser{i}(:,v),f3);
                FR1_nolaser{u,i}(:,v) = interp1(f2_nolaser{i}(:,v),nolaser{i}(:,v),f3);
                
                %Normalize interpolated tuning curves:
                normFR1_laser{u,i}(:,v) = FR1_laser{u,i}(:,v)./max(FR1_nolaser{u,i}(:,v));
                normFR1_nolaser{u,i}(:,v) = FR1_nolaser{u,i}(:,v)./max(FR1_nolaser{u,i}(:,v));
                
            end
        end        

    end
end

a = reshape(normFR1_laser,1,numel(normFR1_laser));
catFR1_laser = [a{:}];

a = reshape(normFR1_nolaser,1,numel(normFR1_nolaser));
catFR1_nolaser = [a{:}];

meanFR1_laser = nanmean( catFR1_laser ,2 );
meanFR1_nolaser = nanmean(catFR1_nolaser,2);

figure; plot(f2,meanFR1_nolaser,'k')
hold on; plot(f2,meanFR1_laser,'r');
hold off; box off;
ylabel('Normalized FR')
xlabel('Octaves from best frequency')

%% Good example tuning curves (so calculate and save images from TC_Select3)
load('C:\Users\Jennifer\Dropbox\GeffenLab\Jennifer\SfN2015\ttestGoodCells2') %Load file names for good cells selected for analysis
for u = 1:length(MNum)
    for i = 1:length(Filters) - 1
        h1 = GOODCELLS{u};
        TC_Select3(MNum(u),h1,Filters(i),1);
    end
end

