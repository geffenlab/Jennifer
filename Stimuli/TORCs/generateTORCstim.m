%% Generate individual TORCs

%Parameters (inputs for makeTorc/makeRipple functions)
omega = [0:0.2:1.4]; % spectral density cycles/octave
w = [4:4:48]; % angular freq Hz
duration = .5; % in seconds
fs = 400000; % sample rate
d = 1; % modulation depth % of mean 
ntones = 501; % number of tones
fmin = 5000; % min freq
oct = 3; % number of octaves

%Create TORCs

nTORC1 = length(omega); %Number of TORCs
for ii = 1:nTORC1
    
    [TORCs{ii},~] = makeTorc(duration,fs,d,omega(ii),w,ntones,fmin,oct);
    
end

w2 = [-4:-4:-48]; % angular freq Hz
for ii = 2:nTORC1
    
    [TORCs{ii + nTORC1 - 1},~] = makeTorc(duration,fs,d,omega(ii),w2,ntones,fmin,oct);
    
end

%% Generate full stimulus 

stimname = 'TORC001';
omegaID = [0:0.2:1.4 -0.2:-0.2:-1.4]; %Negative values denote negative w values, not negative omega values

%Randomly permute order of presenting TORCs
nTORC = length(TORCs);
order = randperm(nTORC);
order_omega = omegaID(order); 
randTORC = TORCs(order);

flp_order = fliplr(order);
flp_order_omega = fliplr(order_omega);
flp_randTORC = TORCs(flp_order);

%Interleave flipped presentation order
tmp_allTORC = [randTORC;flp_randTORC];
tmp_allorder_omega = [order_omega; flp_order_omega];

allTORC = reshape(tmp_allTORC,1, size(tmp_allTORC,1)*size(tmp_allTORC,2));
allorder_omega = reshape(tmp_allorder_omega,1,size(tmp_allorder_omega,1)*size(tmp_allorder_omega,2));

%Add silence in between TORCs
isi = 0.5; %(duration of silence between TORCs in seconds)
isidur = isi*fs;

totaldur = length(allTORC)*duration*fs + length(allTORC)*isi*fs; %Total duration of stimulus

count = 1;
data = zeros(totaldur, 1);
for i = 1:length(allTORC)
    data(count:count+length(allTORC{i})-1) = allTORC{i};
    count = count + length(allTORC{i}) + isidur;
end

data = data./max(abs(data)); %Scale between 0 and 1 for filter to work
data = data./10; %Divide by 10 because current stimulus computer multiplies by 10 (08/16/17)
ttl = zeros(size(data)); %Add second event channel (not used so all can be 0)


% save parameters and stimulus
cd('C:\Users\Jennifer\Documents\MATLAB\Stimuli\TORCs')
save(stimname,'allorder_omega','data','isi','omega','d','w','w2','duration','fs','ntones','fmin','oct','omegaID','TORCs');
wavwrite([data ttl],fs,32,[stimname '.wav']);
