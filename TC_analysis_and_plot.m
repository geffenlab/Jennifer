% Run TC_calculate.m in order to calculate necessary files for this script.
condition = 'CaMKII + ChR2 in IC (anesthetized)';
cd('C:\Users\Jennifer\Documents\MATLAB\Cosyne2018');
allFiles = ls([condition '*.mat']);
AMPS = [3 4;5 6;7 8];
%% Select tone responsive cells with cell quality < 5 and calculate spontaneous and tone-evoked responses

%Pre-allocate variables
spontOFFall = cell(1,size(AMPS,1));
spontONall = cell(1,size(AMPS,1));
magOFFall = cell(1,size(AMPS,1));
magONall = cell(1,size(AMPS,1));
GOODCELL = cell(1,size(AMPS,1));
IDX = cell(1,size(AMPS,1));
baseVARall = cell(1,size(AMPS,1));
OnsetResp1all = cell(1,size(AMPS,1));
OnsetResp2all = cell(1,size(AMPS,1));
OffsetRespall = cell(1,size(AMPS,1));

for n = 1:size(AMPS,1)
    load([condition '_amps ' num2str(AMPS(n,:)) '.mat']);
    
    %Pre-allocate variables
    GOODCELL{n} = cell(1,length(MNum));
    spontIDX = 185:235; %Indices for spontaneous response
    toneIDX = 240:315; %Indices for tone-evoked response
    IDX{n} = cell(1,length(MNum));
    spontON = cell(1,length(MNum));
    spontOFF = cell(1,length(MNum));
    aveON = cell(1,length(MNum));
    aveOFF = cell(1,length(MNum));
    magON = cell(1,length(MNum));
    magOFF = cell(1,length(MNum));
    FRon = cell(1,length(MNum));
    FRoff = cell(1,length(MNum));
    baseVAR = cell(1,length(MNum));
    OnsetResp1 = cell(1,length(MNum));
    OnsetResp2 = cell(1,length(MNum));
    OffsetResp = cell(1,length(MNum));

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
        IDX{n}{u} = find(~isnan(idxGOOD{u}));
        GOODCELL{n}{u} = allCELL{u}(IDX{n}{u},:);

        %Calculate spontaneous firing rate and tone-evoked response magnitude
        numcellB = length(IDX{n}{u});
        spontON{u} = NaN(1,numcellB); spontOFF{u} = NaN(1,numcellB); magON{u} = NaN(1,numcellB);
        magOFF{u} = NaN(1,numcellB); aveON{u} = NaN(1,numcellB); aveOFF{u} = NaN(1,numcellB);
        for j = 1:numcellB   
            spontON{u}(j) = mean(FR_LASER{u}(IDX{n}{u}(j),spontIDX)); spontOFF{u}(j) = mean(FR_NOLASER{u}(IDX{n}{u}(j),spontIDX));
            aveON{u}(j) = mean(FR_LASER{u}(IDX{n}{u}(j),toneIDX)); aveOFF{u}(j) = mean(FR_NOLASER{u}(IDX{n}{u}(j),toneIDX));
            magON{u}(j) = aveON{u}(j) - spontON{u}(j); magOFF{u}(j) = aveOFF{u}(j) - spontOFF{u}(j);
        end

    end

    spontOFFall{n} = [spontOFF{:}];
    spontONall{n} = [spontON{:}];
    magOFFall{n} = [magOFF{:}];
    magONall{n} = [magON{:}];

    NumberOfCells = length([IDX{n}{:}])
    
    % Calculate baseline activity and activity in 3 other time periods for
    % use in next section
    for u = 1:length(MNum)
        if ~isempty(IDX{n}{u})
            FRon{u} = FR_LASER{u}(IDX{n}{u},:);
            FRoff{u} = FR_NOLASER{u}(IDX{n}{u},:);
            baseVAR{u} = std(FRoff{u}(:,spontIDX),0,2); %baseline variability
            OnsetResp1{u} = mean(FRon{u}(:,140:180),2);
            OnsetResp2{u} = mean(FRon{u}(:,200:240),2);
            OffsetResp{u} = mean(FRon{u}(:,390:430),2);
        end
    end

    baseVARall{n} = vertcat(baseVAR{:})';
    OnsetResp1all{n} = vertcat(OnsetResp1{:})';
    OnsetResp2all{n} = vertcat(OnsetResp2{:})';
    OffsetRespall{n} = vertcat(OffsetResp{:})';
    
end

%% Remove cells with laser onset/offset artifacts

Thresh = 3; %Threshold
for n = 1:size(AMPS,1)
    nCell = length(baseVARall{n});
    GOODCELLall{n} = strvcat(GOODCELL{n}{:});
    
    %Categorize cells
    CellCat = NaN(1,nCell);
    for i = 1:nCell
        if OnsetResp2all{n}(i) > spontOFFall{n}(i) + Thresh*baseVARall{n}(i) || OnsetResp2all{n}(i) < spontOFFall{n}(i) - Thresh*baseVARall{n}(i) %Must have significant difference 25 ms before tone onset
            CellCat(i) = 1; % SUSTAINED LASER RESPONSE
        elseif OnsetResp1all{n}(i) > spontOFFall{n}(i) + Thresh*baseVARall{n}(i) || OnsetResp1all{n}(i) < spontOFFall{n}(i) - Thresh*baseVARall{n}(i)
            CellCat(i) = 2; % ONSET RESPONSE 
        elseif OffsetRespall{n}(i) > spontOFFall{n}(i) + Thresh*baseVARall{n}(i) || OffsetRespall{n}(i) < spontOFFall{n}(i) - Thresh*baseVARall{n}(i)
            CellCat(i) = 3; % OFFSET ONLY RESPONSE 
        elseif spontONall{n}(i) < spontOFFall{n}(i) + Thresh*baseVARall{n}(i) && spontONall{n}(i) > spontOFFall{n}(i) - Thresh*baseVARall{n}(i)
            CellCat(i) = 0; % NO LASER RESPONSE 
        end    
    end

    %Select only units that are category 0 and category 1
    badIDX = find(CellCat == 2 | CellCat == 3 | isnan(CellCat));
    GOODCELLall{n}(badIDX,:) = [];

    spontOFFall2{n} = spontOFFall{n};
    spontOFFall2{n}(badIDX) = [];
    spontONall2{n} = spontONall{n};
    spontONall2{n}(badIDX) = [];

    magOFFall2{n} = magOFFall{n};
    magOFFall2{n}(badIDX) = [];
    magONall2{n} = magONall{n};
    magONall2{n}(badIDX) = [];

end

%% Calculate sparseness

nFreq = 50;
SidxON = cell(1,size(AMPS,1));
SidxOFF = cell(1,size(AMPS,1));

for n = 1:size(AMPS,1)
    numcell = size(GOODCELLall{n},1);
    FR_LASER_ALL = [];
    FR_NOLASER_ALL = [];

    h1 = GOODCELLall{n};
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

    %Calculate and smooth firing rate (NOTE: difference from TC_calculate
    % is that this includes all frequencies and all amplitudes)
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
    SidxON{n} = nan(1,numcell);
    SidxOFF{n} = nan(1,numcell);
    for v = 1:numcell
       for j = 1:nFreq
           mFR_ON(j) = mean(FR_LASER_ALL(j,v,241:290));
           mFR_OFF(j) = mean(FR_NOLASER_ALL(j,v,241:290));
       end

        SidxON{n}(v) = Sparseness(mFR_ON,nFreq);  
        SidxOFF{n}(v) = Sparseness(mFR_OFF,nFreq);  
    end

end

%% Produce plots
fontname = 'Arial';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'DefaultTextColor','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spontaneous and tone-evoked response magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for n = 1:size(AMPS,1)
    subplot(size(AMPS,1),4,1+4*(n-1)); scatter(spontOFFall2{n}, spontONall2{n},25,'filled')
    LIM = max([spontOFFall2{n}, spontONall2{n}]);
    %LIM = 40;
    hold on; line([0 LIM], [0 LIM],'Color','k','LineStyle','--');
    hold off;
    set(gca,'TickDir','out'); axis square
    set(gca,'xlim',[0 LIM],'ylim',[0 LIM])
    xlabel('LASER OFF')
    ylabel('LASER ON')
    title('Spontaneous activity (Hz)')

    subplot(size(AMPS,1),4,2+4*(n-1)); bar([mean(spontOFFall2{n}) mean(spontONall2{n})],0.5,'EdgeColor','none');
    hold on; errorbar([mean(spontOFFall2{n}) mean(spontONall2{n})],[std(spontOFFall2{n})./sqrt(length(spontOFFall2{n})) std(spontONall2{n})./sqrt(length(spontONall2{n}))],'k','LineStyle','none')
    box off
    ylabel('Firing rate (Hz)')
    set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
    %set(gca,'YLim',[0 12])
    [p_spont, h_spont] = signrank(spontOFFall2{n}, spontONall2{n});

    subplot(size(AMPS,1),4,3+4*(n-1)); scatter(magOFFall2{n}, magONall2{n},25,'filled')
    LIM = max([magOFFall2{n}, magONall2{n}]);
    %LIM = 100;
    hold on; line([0 LIM], [0 LIM], 'Color','k','LineStyle','--');
    hold off;
    set(gca,'TickDir','out'); axis square
    set(gca,'xlim',[0 LIM],'ylim',[-20 LIM])
    xlabel('LASER OFF')
    ylabel('LASER ON')
    title('Tone-evoked response magnitude (Hz)')

    subplot(size(AMPS,1),4,4+4*(n-1)); bar([mean(magOFFall2{n}) mean(magONall2{n})],0.5,'EdgeColor','none');
    hold on; errorbar([mean(magOFFall2{n}) mean(magONall2{n})],[std(magOFFall2{n})./sqrt(length(magOFFall2{n})) std(magONall2{n})./sqrt(length(magONall2{n}))],'k','LineStyle','none')
    box off
    ylabel('Firing rate (Hz)')
    set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out'); axis square
    %set(gca,'YLim',[0 30])
    [p_mag, h_mag] = signrank(magOFFall2{n}, magONall2{n});
    
end

suptitle(condition);

%%%%%%%%%%%%%
% Sparseness
%%%%%%%%%%%%%

figure; 
for n = 1:size(AMPS,1)
    subplot(size(AMPS,1),2,1+2*(n-1));
    scatter(SidxOFF{n}, SidxON{n}, 25,'filled')
    hold on; line([0 1], [0 1],'Color','k','LineStyle','--')
    hold off; axis square
    xlabel('Sparseness OFF')
    ylabel('Sparseness ON')
    title('Sparseness')
    set(gca,'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1])
    [p_sparse, h_sparse] = signrank(SidxOFF{n}, SidxON{n})

    subplot(size(AMPS,1),2,2+2*(n-1));
    bar([nanmean(SidxOFF{n}) nanmean(SidxON{n})],0.5,'EdgeColor','none');
    hold on; errorbar([nanmean(SidxOFF{n}) nanmean(SidxON{n})],[nanstd(SidxOFF{n})./sqrt(length(SidxOFF{n})) nanstd(SidxON{n})./sqrt(length(SidxON{n}))],'k','LineStyle','none')
    box off
    ylabel('Sparseness')
    set(gca, 'XTickLabels',{'OFF';'ON'},'TickDir','out','YLim',[0 1]); axis square
end

suptitle(condition)