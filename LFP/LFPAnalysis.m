%% Read in mouse/sessions of interest with stimulus indices from ETS file
opsin = 'chr2';
mouseline = 'pv';
loc = 'ac';
cond = 'all';

[MNum, Sesh, TITLE, TUNEidx] = loadsessions(opsin,mouseline,loc,cond);



%%

%%%%%%%%%%%%%%%%%%%%%%% 0. Load and downsample LFP data %%%%%%%%%%%%%%%%%%%%%%%
for mm = 6%3:length(MNum)
    for ss = 2%1:length(Sesh{mm})
        cd(['D:\Spikes\M' num2str(MNum(mm)) '\Session' num2str(Sesh{mm}(ss))]);
        mkdir('LFPAnalysis');

        [TSEv, EventIDs, Nttls, Extras, EventStrings] = Nlx2MatEV('data\Events.nev', [1 1 1 1 1], 0, 1);
        startind= TUNEidx{mm}(ss,1);
        endind = TUNEidx{mm}(ss,2);
        recTS = [TSEv(startind) TSEv(endind)];

        %Index timestamps for
        [ts,header] = Nlx2MatCSC('data\CSC1.ncs',[1 0 0 0 0],1,1);
        startind = find(ts<=recTS(1),1,'last');
        endind = find(ts>=recTS(2),1,'first');
        recind = [startind endind]-1; % zero indexed: tried on 7/24/18 to address offset
        ns = length(recind(1):recind(2));

        chunkSize = 1e4;
        chunkSamps = chunkSize * 512;
        nChunks = ceil(ns / chunkSize);

        nChan = 32;
        RawTS = [];
        fsSample = 2000;
        win = 0.1; %Time around tone onset to extract data (in seconds)
        load(['TS_M' num2str(MNum(mm)) '_Session' num2str(Sesh{mm}(ss)) '.mat']); %Contains Events variable calculated from ETS file
        evlen = cellfun(@length,Events{1});
        EventN = find(evlen == 400); %Event block numbers of interest in Events (mine are all 400 in length)
        toneLFP = cell(1,nChan);
        for i = 1:nChan
            fprintf('Loading channel %d \n', i);
            fn = ['data\CSC' num2str(i) '.ncs'];
            data = [];
            % for each chunk
            for t = 1:nChunks
                if t == nChunks
                    % for last chunk
                    ind = [(recind(1)+chunkSize*(t-1)) recind(2)];
                else
                    %ind = [(chunkSize*(t-1))+recind(1) chunkSize*t];
                    ind = [(chunkSize*(t-1)) chunkSize*t-1] +recind(1);
                end
                fprintf('\tChunk %d/%d -- index = [%d %d]: \n',t,nChunks,ind(1),ind(2));

                tic;
                fprintf('reading... ');
                % preallocate data matrix
                FieldSelection(1) = 1;%timestamps
                FieldSelection(2) = 0;
                FieldSelection(3) = 1;%sample freqs
                FieldSelection(4) = 0;
                FieldSelection(5) = 1;

                ExtractHeader = 0;
                ExtractMode = 2;

                %Load data and timestamps
                if i == 1 %For first channel also pull out timestamps and samping rate
                    [ts, samplefreq, samples] = Nlx2MatCSC(fn, FieldSelection, ExtractHeader, ExtractMode, ind);        
                    data = [data reshape(samples,[],512*(ind(2)-ind(1)+1))];
                    RawTS = [RawTS ts];
                    fsRaw = samplefreq(1);

                    %Downsample timestamps
                    interpdiff = mean(diff(RawTS))/512;
                    interpTSmat = [];
                    for t = 1:length(RawTS)
                        a = repmat(interpdiff,1,511);
                        b = 1:511;
                        c = repmat(RawTS(t),1,511);
                        stepTS = a.*b + c;

                        interpTSmat(:,t) = [RawTS(t); stepTS'];
                    end
                    dur = length(data)/fsRaw;
                    interpTS = reshape(interpTSmat,[],512*length(RawTS));
                    SampleTS = interp1(0:1/fsRaw:(dur-1/fsRaw),interpTS,0:1/fsSample:(dur - 1/fsSample));

                else %For remaining only need to extract data
                    FieldSelection(1) = 0;
                    FieldSelection(3) = 0;
                    samples = Nlx2MatCSC(fn, FieldSelection, ExtractHeader, ExtractMode, ind);
                    data = [data reshape(samples,[],512*(ind(2)-ind(1)+1))];
                end

            end
            RawLFP = data;
            dur = length(RawLFP)/fsRaw;

            fprintf('downsampling data... \n')
            SampleLFP = interp1(0:1/fsRaw:(dur-1/fsRaw),RawLFP,0:1/fsSample:(dur - 1/fsSample));

            %%%%%%%%%%%%%%%%%%%%%%% 1. Pull out relevant data %%%%%%%%%%%%%%%%%%%%%%%  
            %Doing inside the loop because Gogol ran out of memory saving raw or
            %downsampled data

            fprintf('extracting data... \n')
            for u = EventN
                for v = 2:length(Events{1}{u})
                    [~, matchTS] = min(abs(SampleTS - Events{1}{u}(v))); %match Event timestamps with LFP timestamps (Laser OFF trials only)
                    toneLFP{i} = [toneLFP{i}; SampleLFP((matchTS - fsSample*win):(matchTS + fsSample*win))]; %Extract window of data around tone onset
                end
            end

            clear SampleLFP RawLFP
            data = [];
        end
        save('LFPAnalysis\CSDAnalysis','toneLFP','fsSample','EventN');

        %%%%%%%%%%%%%%%%%%%%%%% 2. Average LFPs %%%%%%%%%%%%%%%%%%%%%%%
        % Average across trials, then across triodes, so have a mean LFP per triode

        for i = 1:nChan
           meanLFP(i,:) = mean(toneLFP{i},1); 
        end

        nTT = 10; %Number of triode clusters
        for i = 1:nTT
            ind = (i-1)*3 + 2; %Drop first and last channel
            TTmeanLFP(i,:) = (-1)*mean(meanLFP(ind:ind+2,:),1);
        end

        %%%%%%%%%%%%%%%%%%%%%%% 3. Compute CSD %%%%%%%%%%%%%%%%%%%%%%%
        TT_CSD = [];
        for t = 1:size(TTmeanLFP,2)
            TT_CSD(:,t) = conv(TTmeanLFP(:,t),[1 -2 1],'valid');
        end

        save('LFPAnalysis\CSDAnalysis','TTmeanLFP','TT_CSD','-append');

    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%% 4. Plot CSD %%%%%%%%%%%%%%%%%%%%%%%

for mm = 1:length(MNum)
    for ss = 1:length(Sesh{mm})
        cd(['D:\Spikes\M' num2str(MNum(mm)) '\Session' num2str(Sesh{mm}(ss)) '\LFPAnalysis']);
        load('CSDAnalysis')
        plotCSD(TT_CSD,0,-200:200,2:9);
        suptitle(['M' num2str(MNum(mm)) ' Session ' num2str(Sesh{mm}(ss))]);
        print(['CSD_M' num2str(MNum(mm)) '_' num2str(Sesh{mm}(ss))],'-r400','-djpeg');
    end
end

