%% ************************************************************************
%  *****                   1. EXPERIMENT SESSIONS                     *****
%  ************************************************************************
clear all
opsin = 'chr2';
loc = 'ic';

[MNum, ~, ~] = loadsessions2(opsin,loc);

%% ************************************************************************
%  *****                 2. MAKE PLOTS FOR TUNING                     *****
%  ************************************************************************

%Stimulus parameters
tLaserON = [0.4 0.48 0.508];
tLaserOFF = [0.65 0.73 0.758];

load('D:\Code\TuningCurve\TC003_170815LEFT_stimInfo');
fList = stimInfo.index(:,1);
aList = stimInfo.attenuations;
        
amps = 1;
nAmps = length(amps);
StimOrder = [fList(stimInfo.order)'; aList(ones(length(stimInfo.order),1))'];
trialdur = (stimInfo.tDur+stimInfo.ITI)/1000;
dur = trialdur*length(stimInfo.order);
Time = 0:trialdur:dur; %Each repetition is 400 seconds long, each trial is 500ms long
Time = Time(1:end-1);
StimOrder_Laser = [Time(2:2:end); StimOrder(:,2:2:end)];
StimOrder_NoLaser = [Time(1:2:end); StimOrder(:,1:2:end)];

Win = 0.24;


figure;
for n = 1%:length(MNum)
    cd(['D:\Spikes\M' num2str(MNum(n)) '\TCs']);
    h = dir('data\TC003*01.mat');
    
    for v = 1%:length(h)
        
% I. Pure tone response timecourse
        for i = 1:length(tLaserON) %Number of different laser conditions
            if exist(['data\' h(v).name(1:end-5) num2str(i) '.mat'],'file')
                load(['data\' h(v).name(1:end-5) num2str(i)]);
                subplot(3, 4, 4*(i-1)+1);
                tt1 = find(TCon.t>=.2, 1, 'first');
                tt5 = find(TCon.t>=.8, 1, 'first');
                fr1 = SmoothGaus(mean(TCon.FRmat, 1), 3);
                fr2 = SmoothGaus(mean(TCoff.FRmat, 1), 3);
                plot(TCoff.t(tt1:tt5), fr1(tt1:tt5), 'k');
                hold on;
                plot(TCoff.t(tt1:tt5), fr2(tt1:tt5), 'r');
                ylabel('Firing rate (Hz)');
                line([tLaserON(i) tLaserOFF(i); tLaserON(i) tLaserOFF(i)], [0 0; max(fr1) max(fr1)], 'color', 'g');
                line([0.5 0.55; 0.5 0.55], [0 0; max(fr1) max(fr1)], 'color', 'k');
                box off; axis tight;
                hold off;
            end
        end
        xlabel('Time (s)')

% II. Time vs frequency
        binsize = 0.001;
        for i = 1:length(tLaserON)
            load(['D:\Spikes\M' num2str(MNum(n)) '\SpikeMat\TC003_LEFT-0' num2str(i) h(v).name(6:end-13) '.mat']); 
            SpkTime_NoLaser = SpikeTime(StimOrder_NoLaser,SpikeData,nRep,Time, Win);
            SpkTime_Laser = SpikeTime(StimOrder_Laser,SpikeData,nRep,Time, Win);

            SpikeDataLaser = cell(1,length(fList));        
            SpikeDataNoLaser = cell(1, length(fList));

            for a = 1:length(fList)
                for b = 1:nRep
                    LASER_SPIKE2 = [sort(SpkTime_Laser{a,1,b}); b*ones(1,length(SpkTime_Laser{a,1,b}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                    SpikeDataLaser{a} = [SpikeDataLaser{a} LASER_SPIKE2];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
                end
                FR_LASER(a,:) = smoothFRx4(SpikeDataLaser{a},nRep*stimInfo.repeats,binsize,[-Win Win],5);
            end
            
            for a = 1:length(fList)
                for b = 1:nRep
                    NOLASER_SPIKE2 = [sort(SpkTime_NoLaser{a,1,b}); b*ones(1,length(SpkTime_NoLaser{a,1,b}))]; %Add a "trial #" so final format will be similar to SpikeData arrays
                    SpikeDataNoLaser{a} = [SpikeDataNoLaser{a} NOLASER_SPIKE2];%Concatenate "trials" to format like SpikeData(3,:) and SpikeData(4,:)
                end
                FR_NOLASER(a,:) = smoothFRx4(SpikeDataNoLaser{a},nRep*stimInfo.repeats,binsize,[-Win Win],5);
            end
            

            PUSH(1) = 0;
            subplot(3,4,[i+1 i+5]); hold on;
            for a = 1:length(fList)
                if a > 1
                    PUSH(a) = max([FR_LASER(a-1,:) FR_NOLASER(a-1,:)]) + PUSH(a-1) + 10;
                end
                plot(FR_NOLASER(a,:) + PUSH(a),'k');
                plot(FR_LASER(a,:) + PUSH(a),'r');
            end
            line([240 240], [0 PUSH(end)],'Color','k','LineStyle','--'); %Tone onset
            line([290 290], [0 PUSH(end)],'Color','k','LineStyle','--'); %Tone offset
            Lon = (tLaserON(i) - 0.5)*1000 + 240;
            Loff = Lon + 250;
            line([Lon Lon], [0 PUSH(end)], 'Color','g') %Laser onset 
            line([Loff Loff], [0 PUSH(end)], 'Color','g') %Laser offset 
            set(gca,'TickDir','out', 'XTick',40:100:480,'XTickLabel',-0.2:0.1:0.2,'YTickLabel',[])
            xlabel('Time  (s)'); ylabel('Frequency')
            hold off;
        end
        
% III. Tuning curve
        for i = 1:length(tLaserON)
            load(['data\' h(v).name(1:end-5) num2str(i)]);
            for u = 1:nAmps
                Ton(:, u) = SmoothGaus(TCon.TCmat{1}(:, amps(u)), 3);
            end

            for u = 1:nAmps
                Toff(:, u) = SmoothGaus(TCoff.TCmat{1}(:, amps(u)), 3);
            end

            subplot(3,4,4*2 + i + 1)
            plot(TCoff.F, mean(Ton(:, amps),2), 'r');
            hold on; plot(TCoff.F, mean(Toff(:, amps),2), 'k');
            set(gca, 'xscale', 'log'); 
            xlabel('Frequency (Hz)')
            ylabel('Firing rate (Hz)')
            box off; axis tight;
            %legend('laser', 'no laser','location', 'bestoutside');
            hold off;
        end


        suptitle(['Frequency response: ' h(v).name(1:end-13) ' ' num2str(CellInfo(end))])
        set(gcf,'PaperPositionMode','auto');         
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperUnits','points');
        set(gcf,'PaperSize',[2000 1000]);
        set(gcf,'Position',[0 0 2000 1000]);
%        print('-djpeg',['pics\Tuning-' h(v).name(7:end-13)]);

    end
end

%% ************************************************************************
%  *****                  3. MAKE PLOTS FOR DRCs                      *****
%  ************************************************************************

LaserOn1 = 0.5;
LaserDur = 0.25;
stimdur = 300;

win0 = [0 0.25];
win1 = [0.25 0.5];
win2 = [0.5 0.75];
win3 = [0.75 1];

nFiles = 8;
StimParams = cell(1,8);
for i = 1:nFiles
    StimParams{i} = load(['DRC001-0' num2str(i)],'params'); %Load stimulus parameters
end

for n = 1:length(MNum)
    cd(['D:\Spikes\M' num2str(MNum(n))]);
    h = dir('SpikeMat\DRC001_LEFT-01*.mat');

    nFiles = 8;

    for v = 1:length(h)
        for i = 1:nFiles

            load(['SpikeMat\DRC001_LEFT-0' num2str(i) '-' h(v).name(16:end)]); %Load spike times 

            spikes = SpikeData(3,:);
            spikes(spikes >= stimdur) = [];
            SpikeTmod = mod(spikes,1);
            SpikeFR(i,:) = smoothFRx4(SpikeTmod,stimdur,0.001,[0 1],5);
        end

        subplot(3,2,[1 2]);
        plot(mean(SpikeFR,1),'k')
        line([LaserOn1*1000 LaserOn1*1000],[0 max(mean(SpikeFR,1))],'Color','g')
        line([(LaserOn1+LaserDur)*1000 (LaserOn1+LaserDur)*1000], [0 max(mean(SpikeFR,1))], 'Color','g');
        set(gca,'TickDir','out','XTick',0:100:1000,'XTickLabel',0:0.1:1)
        xlabel('Time (s)'); ylabel('Firing rate (Hz)')
        box off; axis tight;


        load(['TCs\data\DRC001-' h(v).name(16:end)])
        %Plot STRFs
        MM = max(max([STAon;STAoff1;STAoff2;STAoff3]));
        mm = min(min([STAon;STAoff1;STAoff2;STAoff3]));

        suptitle(['DRCs (DRC001): ' h(v).name(16:end-4) ' ' num2str(CellInfo(end))])
        subplot(3,2,3); plotSTA([-0.1:0.005:0],StimParams{1}.params.freqs/1000,STAon,1,[mm,MM]); colorbar;
        title('Laser On')
        subplot(3,2,4); plotSTA([-0.1:0.005:0],StimParams{1}.params.freqs/1000,STAoff1,1,[mm,MM]); colorbar;
        title('Laser Off 0.25 - 0.5')
        subplot(3,2,5); plotSTA([-0.1:0.005:0],StimParams{1}.params.freqs/1000,STAoff2,1,[mm,MM]); colorbar;
        title('Laser Off 0.5 - 0.75')
        subplot(3,2,6); plotSTA([-0.1:0.005:0],StimParams{1}.params.freqs/1000,STAoff3,1,[mm,MM]); colorbar;
        title('Laser Off 0.75 - 1')

        set(gcf,'PaperPositionMode','auto');         
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperUnits','points');
        set(gcf,'PaperSize',[1000 1000]);
        set(gcf,'Position',[0 0 1000 1000]);

        print('-djpeg',['TCs\pics\DRC-'  h(v).name(16:end-4) '_STRF.jpg']);
    end
end