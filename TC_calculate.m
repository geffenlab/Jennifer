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
TITLE = 'CaMKII + ChR2 in IC (anesthetized)';
MNum = [3128 3129 3130 3140 3141];
Sesh = cell(1,length(MNum));
Sesh{1} = [2,3,4]; Sesh{2} = [2,3,4]; Sesh{3} = [2,3,4]; Sesh{4} = [2,3];
Sesh{5} = [2,3];
 
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

%% Average tuning curve stimulus across top two sound levels

nfreq = 7;
amps = [7,8];

FR_NOLASER = cell(1,length(MNum));
FR_LASER = cell(1,length(MNum));
CellQ = cell(1,length(MNum));
for u = 1:length(MNum)  
    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    h1 = ls('data\TC3-1*_laser.mat');
    idxUSE = find(ismember(str2num(h1(:,15)),Sesh{u}));
    h = h1(idxUSE,:);
    for v = 1:size(h,1)
        [FR_NOLASER{u}(v,:), FR_LASER{u}(v,:), CellQ{u}(v)] = smoothTC(MNum(u),h(v,:),nfreq,amps);
    end
    
    allCELL{u} = h;
end

cd('C:\Users\Jennifer\Documents\MATLAB\Cosyne2018');
save([TITLE '_amps ' num2str(amps)],'TITLE','MNum','Sesh','allCELL','FR_LASER','FR_NOLASER','CellQ','nfreq','amps');

%% Average tuning curve stimulus across 40/50 dB

nfreq = 7;
amps = [5,6];

FR_NOLASER = cell(1,length(MNum));
FR_LASER = cell(1,length(MNum));
CellQ = cell(1,length(MNum));
for u = 1:length(MNum)  
    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    h1 = ls('data\TC3-1*_laser.mat');
    idxUSE = find(ismember(str2num(h1(:,15)),Sesh{u}));
    h = h1(idxUSE,:);
    for v = 1:size(h,1)
        [FR_NOLASER{u}(v,:), FR_LASER{u}(v,:), CellQ{u}(v)] = smoothTC(MNum(u),h(v,:),nfreq,amps);
    end
    
    allCELL{u} = h;
end

cd('C:\Users\Jennifer\Documents\MATLAB\Cosyne2018');
save([TITLE '_amps ' num2str(amps)],'TITLE','MNum','Sesh','allCELL','FR_LASER','FR_NOLASER','CellQ','nfreq','amps');

%% Average tuning curve stimulus across 20/30 dB

nfreq = 7;
amps = [3,4];

FR_NOLASER = cell(1,length(MNum));
FR_LASER = cell(1,length(MNum));
CellQ = cell(1,length(MNum));
for u = 1:length(MNum)  
    cd(['D:\Spikes\M' num2str(MNum(u)) '\TCs']);
    h1 = ls('data\TC3-1*_laser.mat');
    idxUSE = find(ismember(str2num(h1(:,15)),Sesh{u}));
    h = h1(idxUSE,:);
    for v = 1:size(h,1)
        [FR_NOLASER{u}(v,:), FR_LASER{u}(v,:), CellQ{u}(v)] = smoothTC(MNum(u),h(v,:),nfreq,amps);
    end
    
    allCELL{u} = h;
end

cd('C:\Users\Jennifer\Documents\MATLAB\Cosyne2018');
save([TITLE '_amps ' num2str(amps)],'TITLE','MNum','Sesh','allCELL','FR_LASER','FR_NOLASER','CellQ','nfreq','amps');