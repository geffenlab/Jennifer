function [output] = FMsweepTC2(FileName, win, binsize, smoothWin)
%Calculate tuning of FM sweep stimulus J002
%
%Inputs:
%   FileName = string with full file name of SpikeData
%   win = Where (relative to tone onset) to start looking for spikes (in seconds).
%   binsize = size of time bins for smoothing (in seconds)
%   smoothWin = Number of bins over which to smooth (using SmoothGaus
%               function)
%
%Outputs:
%   tuneIDX = FM tuning index (Sparseness) [IDXoff IDXon]
%   dirIDX = FM sweep direction selectivity index (see Carruthers et al
%            2013) [IDXoff IDXon]

load('D:\Code\TuningCurve\J002F105') % Load FM sweep stimulus (although labeled F105, this does not change stimulus order so it doesn't matter F105 or F104)
load(FileName);  

t = [0:isi:length(x)/fs]; %Start times for each sweep
t = t(1:end-1);

if ~isempty(SpikeData)
    %Pull out spikes for each trial (FULL TRIAL) across all repetitions
    ResponseT = cell(nRep,length(t));
    SpikeT = SpikeData(2,:);

    for i = 1:nRep
        %a = find(SpikeData(4,:) == i); %Find spikes that occur in repetition i
        %SpikeT = SpikeData(3,a); %Time of spikes that occur within repetition i                   
        for j = 1:length(t)
            %b = find(SpikeT < t(j) + isi & SpikeT > t(j) ); 
            b = find(SpikeT < t(j) +win + isi +length(t)*(i-1) & SpikeT > t(j) +win + length(t)*(i-1) ); 
            if ~isempty(b)
                ResponseT{i,j} = SpikeT(b) - t(j) - length(t)*(i-1); %Time of spike occurance after tone onset.
            end
        end   
    end 

    %Allocate the first stim of first rep which cannot look before onset...
    c = find(SpikeT > t(end) +win + isi + length(t)*nRep & SpikeT < t(end) + isi + length(t)*nRep); %Take from end of last stim of last rep
    d = find(SpikeT > t(1) & SpikeT < t(1) + win + isi); 
    if ~isempty(c) 
        tmp_c = SpikeT(c) - t(end) - isi;
    else tmp_c = [];
    end
    ResponseT{1,1} = [tmp_c SpikeT(d)];

    Stims = unique(opersorder); %List of sweep speeds used in stimulus
    StimsLen = NaN(1,length(Stims)); %List of sweep speed lengths (in seconds)
    for i = 1:length(Stims)
        z = find(opersorder == Stims(i));
        StimsLen(i) = lenorder(z(1));
    end


    %Separate LASER ON/LASER OFF trials for each sweep
    ResponseT_NoLaser = ResponseT(:,1:2:end);
    ResponseT_Laser = ResponseT(:,2:2:end);

    StimOrder_NoLaser = opersorder(1:2:end);
    StimOrder_Laser = opersorder(2:2:end);

    SpkTime_NoLaser = cell(length(Stims),nRep);
    for i = 1:length(StimOrder_Laser)
        a = find(Stims == StimOrder_NoLaser(i));
        for j = 1:nRep
           SpkTime_NoLaser{a,j} = ResponseT_NoLaser{j,i};
        end 
    end

    SpkTime_Laser = cell(length(Stims),nRep);
    for i = 1:length(StimOrder_Laser)
        a = find(Stims == StimOrder_Laser(i));
        for j = 1:nRep    
            SpkTime_Laser{a,j} = ResponseT_Laser{j,i};
        end
    end

    % Calculate smoothed FR for each stimulus
    SpikeDataNoLaser = cell(1,length(Stims));
    SpikeDataLaser = cell(1,length(Stims));
    FR_NOLASER = []; FR_LASER = [];
    for i = 1:length(Stims)
        for j = 1:nRep
            tmp1 = [SpkTime_NoLaser{i,j}; j*ones(1,length(SpkTime_NoLaser{i,j}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
            tmp2 = [SpkTime_Laser{i,j}; j*ones(1,length(SpkTime_Laser{i,j}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
            SpikeDataNoLaser{i} = [SpikeDataNoLaser{i} tmp1];           
            SpikeDataLaser{i} = [SpikeDataLaser{i} tmp2];
        end
        FR_NOLASER(i,:) = smoothFRx4(SpikeDataNoLaser{i},nRep,binsize,[win win+isi],smoothWin);
        FR_LASER(i,:) = smoothFRx4(SpikeDataLaser{i},nRep,binsize,[win win+isi],smoothWin);
    end

    %PLOT FIRING RATE FOR EACH STIMULUS FOR FULL TRIAL
    LIM = max([FR_NOLASER(:); FR_LASER(:)]); %Find max FR to scale images
    offset = abs(win(1)/binsize); %Bin offset for plotting stimulus locations

        %figure; 
        subplot(3,2,1);
        imagesc(FR_NOLASER, [0 LIM+1]);
        hold on 
        for i = 1:length(Stims)
            [~, StimBin(i)] = min(abs([0:binsize:isi] - StimsLen(i))); %Length of stimulus (in bins)
            line([StimBin(i)+offset StimBin(i)+offset], [i-0.5 i+0.5],'Color','k'); 
        end
        line([offset offset], [0.5 length(Stims)+0.5],'Color','k'); %Mark tone onset
        hold off        
        set(gca,'TickDir','out', 'YDir','normal','YTick', 2:2:28,'YTickLabel',[-7.5:1:-1.5 2:1:8],'XTick',0.1/binsize:0.1/binsize:isi/binsize, ...
            'XTickLabel',-0.1:0.1:1.9)
        title('No laser'); xlabel('Time from sweep onset (s)'); ylabel('FM sweep speed')
        %colorbar

        %figure; 
        subplot(3,2,2);
        imagesc(FR_LASER, [0 LIM+1]);
        hold on 
        for i = 1:length(Stims)
            line([StimBin(i)+offset StimBin(i)+offset], [i-0.5 i+0.5],'Color','k'); 
        end
        line([offset offset], [0.5 length(Stims)+0.5],'Color','k');
        line([offset - 0.1/binsize offset - 0.1/binsize], [0.5 length(Stims)+0.5], 'Color','r') %Laser onset (100ms prior to sweep onset)
        line([offset + 1.4/binsize offset + 1.4/binsize], [0.5 length(Stims)+0.5], 'Color','r') %Laser offset (1.5s duration)
        hold off        
        set(gca,'TickDir','out', 'YDir','normal','YTick', 2:2:28,'YTickLabel',[-7.5:1:-1.5 2:1:8],'XTick',0.1/binsize:0.1/binsize:isi/binsize,'XTickLabel',-0.1:0.1:1.9)
        title('Laser'); xlabel('Time from sweep onset (s)'); ylabel('FM sweep speed')
        colorbar

        %figure;
        subplot(3,2,[3 5])
        hold on;
        PUSH(1) = 0;
        for i = 1:length(Stims)
            if i > 1
                PUSH(i) = max([FR_LASER(i-1,:) FR_NOLASER(i-1,:)]) + PUSH(i-1) + 10;
            end
            plot(FR_NOLASER(i,:) + PUSH(i),'k');
            plot(FR_LASER(i,:) + PUSH(i),'r');
        end
        line([offset offset], [0 PUSH(end)],'Color','k','LineStyle','--');
        line([offset - 0.1/binsize offset - 0.1/binsize], [0 PUSH(end)], 'Color','r','LineStyle','--') %Laser onset (100ms prior to sweep onset)
        line([offset + 1.4/binsize offset + 1.4/binsize], [0 PUSH(end)], 'Color','r','LineStyle','--') %Laser offset (1.5s duration)
        set(gca,'TickDir','out', 'XTick',0.1/binsize:0.1/binsize:isi/binsize,'XTickLabel',-0.1:0.1:1.9)
        xlabel('Time from sweep onset (s)'); ylabel('Firing Rate (Hz)')

        %Calculate and PLOT max and mean response during sweep presentation +
        %50 ms and subtract NO LASER baseline FR
        PostWin = 0.05; %Window after sweep offset in which to look for spikes (s)
        BaseWin = 0.1; %Window to use as baseline
        BASELINE = FR_NOLASER(:,offset - BaseWin/binsize + 1:offset - 1);
        meanSPONT = mean(BASELINE(:));
        stdSPONT = std(BASELINE(:));
        for i = 1:length(Stims)
            meanOFF(i) = mean(FR_NOLASER(i,offset:StimBin(i) + offset + PostWin/binsize)) - mean(FR_NOLASER(i,offset - BaseWin/binsize:offset - 1));
            meanON(i) = mean(FR_LASER(i,offset:StimBin(i) + offset + PostWin/binsize))- mean(FR_NOLASER(i,offset - BaseWin/binsize:offset - 1));
            devOFF(i) = mean(FR_NOLASER(i,offset:StimBin(i) + offset + PostWin/binsize))./stdSPONT; %Calculate number of standard deviations during sound compared to baseline activity
        end


        %figure; 
        subplot(3,2,4);
        scatter(1:length(Stims),meanOFF,15,'k','filled');
        hold on; plot(SmoothGaus(meanOFF,2),'k'); scatter(1:length(Stims),meanON,15,'r','filled'); plot(SmoothGaus(meanON,2),'r')
        scatter(1:length(Stims),devOFF,15,'b','filled');
        set(gca,'TickDir','out', 'XTick', 2:2:28,'XTickLabel',[-7.5:-1.5 2:8])
        ylabel('Mean FR'); xlabel('FM sweep speed')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate from Isaac's paper
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        binsize2 = 0.01;
        offset2 = abs(win(1)/binsize2);
        t = [win:binsize2:win+isi];
        for i = 1:length(Stims)
            [~, StimBin2(i)] = min(abs([0:binsize2:isi] - StimsLen(i))); %Length of stimulus (in bins)

        end
        %STEP 1: Bin into 10ms bins

        for i = 1:length(Stims)
            f_nolaser = histc(SpikeDataNoLaser{i}(1,:), t)/(nRep*binsize2);
            f_laser = histc(SpikeDataLaser{i}(1,:), t)/(nRep*binsize2);
            f1_nolaser(i,:) = f_nolaser(1:end-1);
            f1_laser(i,:) = f_laser(1:end-1);
        end

        %STEP 2: Take mean and std for 200ms before sweep onset

        f1_spont = f1_nolaser(:,1:offset2-1);
        std_OFF = std(f1_spont(:));
        mean_OFF = mean(f1_spont(:));

        %STEP 3: Find bins for which spike count exceeds 95% confidence limit
        %of gaussian (1.96*STD) and take average and normalize
        ConfLim = 1.96.*std_OFF;
        for i = 1:length(Stims)
            indbins = offset2:StimBin2(i) + PostWin/binsize2 + offset2; %Indices of sweep + 50ms post-sweep
            idx_OFF{i} = find(f1_nolaser(i,indbins) > mean_OFF + ConfLim | f1_nolaser(i,indbins) < mean_OFF - ConfLim);
            idx_OFF{i} = idx_OFF{i} + offset2 - 1;
            idx_ON{i} = find(f1_laser(i,indbins) > mean_OFF + ConfLim | f1_laser(i,indbins) < mean_OFF - ConfLim);
            idx_ON{i} = idx_ON{i} + offset2 - 1;
        end

        normFR_OFF = zeros(1,length(Stims));
        normFR_ON = zeros(1,length(Stims));
        for i = 1:length(Stims)
            if ~isempty(idx_OFF{i}) 
                meanFR_OFF = mean(f1_nolaser(i,idx_OFF{i}));
                normFR_OFF(i) = (meanFR_OFF - mean_OFF);
            end
            if ~isempty(idx_ON{i})
                meanFR_ON = mean(f1_laser(i,idx_ON{i}));
                normFR_ON(i) = (meanFR_ON - mean_OFF);
            end
        end

    if mean(meanOFF) > 1 && max(devOFF) > 2 %Requirement that mean of baseline subtracted FR and at least one sweep have variance of FR over baseline  
    
    % Calculate FM tuning index (sparsensess)
    N_off = length(Stims);
    N_on = length(Stims);
    
    tuneIDX(1) = Sparseness(normFR_OFF,N_off);
    tuneIDX(2) = Sparseness(normFR_ON,N_on);
    
    %Calculate FM directionality index
    
    downOFF = mean(normFR_OFF(1:14));
    upOFF = mean(normFR_OFF(15:end));
    downON = mean(normFR_ON(1:14));
    upON = mean(normFR_ON(15:end));
    
    dirIDX(1) = (upOFF - downOFF)/(upOFF + downOFF);
    dirIDX(2) = (upON - downON)/(upON + upOFF);
    else dirIDX = [NaN NaN]; tuneIDX = [NaN NaN]; tuneIDXboot = [NaN NaN]; dirIDXboot = [NaN NaN];
    end
    
    output.dirIDX = dirIDX;
    output.tuneIDX = tuneIDX;
    output.normFR_ON = normFR_ON;
    output.normFR_OFF = normFR_OFF;
    output.idxON = idx_ON;
    output.idxOFF = idx_OFF;
    output.FR_LASER = FR_LASER;
    output.FR_NOLASER = FR_NOLASER;
    output.meanOFF = meanOFF;
    output.meanON = meanON;
    
    
    subplot(3,2,6);
    scatter(1:length(Stims),normFR_OFF,15,'k','filled');
    hold on; plot(SmoothGaus(normFR_OFF,2),'k'); scatter(1:length(Stims),normFR_ON,15,'r','filled'); plot(SmoothGaus(normFR_ON,2),'r')
    set(gca,'TickDir','out', 'XTick', 2:2:28,'XTickLabel',[-7.5:-1.5 2:8])
    ylabel('norm FR'); xlabel('FM sweep speed');

    suptitle([num2str(CellInfo) ' : tuneOFF = ' num2str(tuneIDX(1)) '; tuneON = ' num2str(tuneIDX(2)) ...
        ' | dirOFF = ' num2str(dirIDX(1)) '; dirON = ' num2str(dirIDX(2))]);
    fig = gcf;
    set(fig,'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'PaperPosition',[0 0 16 10], 'PaperSize',[16 10]);
    
else 
    output.dirIDX = [NaN NaN];
    output.tuneIDX = [NaN NaN];
    output.normFR_ON = NaN;
    output.normFR_OFF = NaN;
    output.idxON = NaN;
    output.idxOFF = NaN;
    output.FR_LASER = NaN;
    output.FR_NOLASER = NaN;
    output.meanOFF = NaN;
    output.meanON = NaN;   
end

