function [MNum,Sesh,TITLE, TUNEidx] = loadsessions(opsin, mouseline, loc, cond);
%Load mouse and session numbers for different experiment conditions
%
%Inputs: (ALL LOWERCASE STRINGS)
%   opsin: name of opsin, 'chr2' or 'archt'
%   mouseline: name of mouseline, 'pv' 'som' or 'camk2'
%   loc: recording location, 'ac' or 'ic
%   cond: anesthesia condition, 'awake' or 'anesthetized'
%
%Outputs:
%   MNum: list of mouse IDs (1xN array)
%   Sesh: list of sessions for each mouse ID (1xN cell array)
%   TITLE: string listing the input experiment info
%   TUNEidx: (optional) contains indices for tuning events, used for LFP
%            analysis
%

if strcmp(cond,'awake')
    if strcmp(opsin,'chr2')
        if strcmp(mouseline,'som')
            if strcmp(loc,'ic')
                TITLE = 'SOM + ChR2 in IC (awake)';
                MNum = [3138 3145 3146 3176 3177 3194 3195];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1,2]; Sesh{2} = [1]; Sesh{3} = [1,2]; Sesh{4} = [1,2]; Sesh{5} = [1,2]; Sesh{6} = [1 2]; Sesh{7} = [1 2];
            elseif strcmp(loc,'ac')
                TITLE = 'SOM + ChR2 in AC (awake)';
                MNum = [3176 3177 3194 3195];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [3]; Sesh{2} = [3,4]; Sesh{3} = [3 4]; Sesh{4} = [3 4];
            end
        elseif strcmp(mouseline,'camk2')
            if strcmp(loc,'ic')                       
                TITLE = 'CaMKII + ChR2 in IC (awake)';
                MNum = [3143 3150 3220 3221];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1]; Sesh{2} = [1,2]; Sesh{3} = [1 2 3 4]; Sesh{4} = [1 2 3 4];
            elseif strcmp(loc,'ac')                
                TITLE = 'CaMKII + ChR2 in AC (awake)';
                MNum = [3220];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [5];
            end
        elseif strcmp(mouseline,'pv')
            if strcmp(loc,'ic')                
                TITLE = 'PV + ChR2 in IC (awake)';
                MNum = [3149 3174 3175 3208];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1,2]; Sesh{2} = [1 2]; Sesh{3} = [1 2]; Sesh{4} = [1 2];
            elseif strcmp(loc,'ac')
                TITLE = 'PV + ChR2 in AC (awake)';
                MNum = [3174 3175];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [3]; Sesh{2} = [3,4];
            end
        end
    elseif strcmp(opsin,'archt') %ArchT
        if strcmp(mouseline,'som')
            if strcmp(loc,'ic')
                TITLE = 'SOM + ArchT in IC (awake)';
                MNum = [3162 3164 3196];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 1; Sesh{2} = [1,2]; Sesh{3} = [1 2];
            elseif strcmp(loc,'ac')
                TITLE = 'SOM + ArchT in AC (awake)';
                MNum = [3161 3162 3196];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1,2]; Sesh{2} = 2; Sesh{3} = [3,4];
            end
        elseif strcmp(mouseline,'camk2')
            if strcmp(loc,'ic')                       
                TITLE = 'CaMKII + ArchT in IC (awake)';
                MNum = [3153 3180 3181 3227 3228 3229];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 2; Sesh{2} = [1,2]; Sesh{3} = [1,2]; Sesh{4} = [1 2 3 4]; Sesh{5} = [1 2 3 4]; Sesh{6} = [1 2 3 4];
            elseif strcmp(loc,'ac')                
                TITLE = 'CaMKII + ArchT in AC (awake)';
                MNum = [3180 3181 3227 3228 3229];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [3,4]; Sesh{2} = 3; Sesh{3} = [5]; Sesh{4} = [5 6]; Sesh{5} = [5 6];
            end
        elseif strcmp(mouseline,'pv')
            if strcmp(loc,'ic')                
                TITLE = 'PV + ArchT in IC (awake)';
                MNum = [3172 3182 3183 3209];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 1; Sesh{2} = [1 2]; Sesh{3} = [1 2]; Sesh{4} = [1 2];
            elseif strcmp(loc,'ac')
                TITLE = 'PV + ArchT in AC (awake)';
                MNum = [3172 3182 3183 3209];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 2; Sesh{2} = [3,4]; Sesh{3} = 3; Sesh{4} = [3 4];

            end
        end
    end
    
    
    
elseif strcmp(cond,'anesthetized')
    if strcmp(opsin,'chr2')
        if strcmp(mouseline,'som')
            if strcmp(loc,'ic')
                TITLE = 'SOM + ChR2 in IC (anesthetized)';
                MNum = [3125 3126 3127 3135 3136 3137];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [2,3]; Sesh{2} = [1,2]; Sesh{3} = [2,3,4]; Sesh{4} = [2 3];
                Sesh{5} = [2 3]; Sesh{6} = [1];
            elseif strcmp(loc,'ac')
                TITLE = 'SOM + ChR2 in AC (anesthetized)';
                MNum = [3127 3135 3136 3158 3159];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1; Sesh{4} = [1,2]; Sesh{5} = 1;
            end
        elseif strcmp(mouseline,'camk2')
            if strcmp(loc,'ic')                       
                TITLE = 'CaMKII + ChR2 in IC (anesthetized)';
                MNum = [3128 3129 3130 3140 3141];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [2,3,4]; Sesh{2} = [2,3,4]; Sesh{3} = [2,3,4]; Sesh{4} = [2,3]; Sesh{5} = [2,3];
            elseif strcmp(loc,'ac')
                TITLE = 'CaMKII + ChR2 in AC (anesthetized)';
                MNum = [3128 3129 3130 3140 3141];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1; Sesh{4} = 1; Sesh{5} = 1;                
            end
        elseif strcmp(mouseline,'pv')
            if strcmp(loc,'ic')                
                TITLE = 'PV + ChR2 in IC (anesthetized)';
                MNum = [3100 3101 3102 3103 3104 3111 3112 3113 3121 3122 3124];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1,2]; Sesh{2} = [1,2]; Sesh{3} = [1,2]; Sesh{4} = [1,2]; Sesh{5} = [1,2];
                Sesh{6} = 1; Sesh{7} = [1,2]; Sesh{8} = 1; Sesh{9} = [2,3]; Sesh{10} = [2,3]; Sesh{11} = 1;
            elseif strcmp(loc,'ac')
                TITLE = 'PV + ChR2 in AC (anesthetized)';
                MNum = [3120 3121 3122 3124];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1]; Sesh{2} = [1]; Sesh{3} = [1]; Sesh{4} = [2];
            end
        end
       
    elseif strcmp(opsin,'archt') %ArchT
        if strcmp(mouseline,'som')
            if strcmp(loc,'ic')
                TITLE = 'SOM + ArchT in IC (anesthetized)';
                MNum = [3139 3144 3163];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [2,3]; Sesh{2} = [2,3]; Sesh{3} = 2;
            elseif strcmp(loc,'ac')
                TITLE = 'SOM + ArchT in AC (anesthetized)';
                MNum = [3139 3144 3160 3163];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1]; Sesh{2} = [1]; Sesh{3} = [1,2]; Sesh{4} = 1;
            end
        elseif strcmp(mouseline,'camk2')
            if strcmp(loc,'ic')                       
                TITLE = 'CaMKII + ArchT in IC (anesthetized)';
                MNum = [3154];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 4;
            elseif strcmp(loc,'ac')                
                TITLE = 'CaMKII + ArchT in AC (anesthetized)';
                MNum = [3154];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [2,3];
            end
        elseif strcmp(mouseline,'pv')
            if strcmp(loc,'ic') 
                TITLE = 'PV + ArchT in IC (anesthetized)';
                MNum = [];
                Sesh = cell(1,length(MNum));
            elseif strcmp(loc,'ac')
                TITLE = 'PV + ArchT in AC (anesthetized)';
                MNum = [];
                Sesh = cell(1,length(MNum));
            end
        end
    end
    
elseif strcmp(cond,'all')
    if strcmp(opsin,'chr2')
        if strcmp(mouseline,'som')
            if strcmp(loc,'ic')
                TITLE = 'SOM + ChR2 in IC';
                MNum = [3125 3126 3127 3135 3136 3137 3138 3145 3146 3176 3177 3194 3195];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [2,3]; Sesh{2} = [1,2]; Sesh{3} = [2,3,4]; Sesh{4} = [2 3];
                Sesh{5} = [2 3]; Sesh{6} = [1]; Sesh{7} = [1,2]; Sesh{8} = [1]; Sesh{9} = [1,2]; 
                Sesh{10} = [1,2]; Sesh{11} = [1,2]; Sesh{12} = [1 2]; Sesh{13} = [1 2];                
            elseif strcmp(loc,'ac')
                TITLE = 'SOM + ChR2 in AC';
                MNum = [3127 3135 3136 3158 3159 3176 3177 3194 3195];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1; Sesh{4} = [1,2]; Sesh{5} = 1;
                Sesh{6} = [3]; Sesh{7} = [3,4]; Sesh{8} = [3 4]; Sesh{9} = [3 4];

            end
        elseif strcmp(mouseline,'camk2')
            if strcmp(loc,'ic')                       
                TITLE = 'CaMKII + ChR2 in IC';
                MNum = [3128 3129 3130 3140 3141 3143 3150 3220 3221];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [2,3,4]; Sesh{2} = [2,3,4]; Sesh{3} = [2,3,4]; Sesh{4} = [2,3]; Sesh{5} = [2,3];
                Sesh{6} = [1]; Sesh{7} = [1,2]; Sesh{8} = [1 2 3 4]; Sesh{9} = [1 2 3 4];
                                
            elseif strcmp(loc,'ac')                
                TITLE = 'CaMKII + ChR2 in AC';
                MNum = [3128 3129 3130 3140 3141 3220];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 1; Sesh{2} = 1; Sesh{3} = 1; Sesh{4} = 1; Sesh{5} = 1; Sesh{6} = [5];
                  
            end
        elseif strcmp(mouseline,'pv')
            if strcmp(loc,'ic')                
                TITLE = 'PV + ChR2 in IC';
                MNum = [3100 3101 3102 3103 3104 3111 3112 3113 3121 3122 3124 3149 3174 3175 3208];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1,2]; Sesh{2} = [1,2]; Sesh{3} = [1,2]; Sesh{4} = [1,2]; Sesh{5} = [1,2];
                Sesh{6} = 1; Sesh{7} = [1,2]; Sesh{8} = 1; Sesh{9} = [2,3]; Sesh{10} = [2,3]; Sesh{11} = 1;
                Sesh{12} = [1,2]; Sesh{13} = [1 2]; Sesh{14} = [1 2]; Sesh{15} = [1 2];

            elseif strcmp(loc,'ac')
                TITLE = 'PV + ChR2 in AC';
                MNum = [3120 3121 3122 3124 3174 3175];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1]; Sesh{2} = [1]; Sesh{3} = [1]; Sesh{4} = [2]; Sesh{5} = [3]; Sesh{6} = [3,4];
                if nargout > 3
                    TUNEidx{1} = [3672 11678]; TUNEidx{2} = [1445 9447]; TUNEidx{3} = [1443 9445];
                    TUNEidx{4} = [11890 19892]; TUNEidx{5} = [1843 9845]; TUNEidx{6} = [1843 9845; 13370 21372];
                end
                
            end
        end
    elseif strcmp(opsin,'archt')
        if strcmp(mouseline,'som')
            if strcmp(loc,'ic')
                TITLE = 'SOM + ArchT in IC';
                MNum = [3139 3144 3163 3162 3164 3196];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [2,3]; Sesh{2} = [2,3]; Sesh{3} = 2; Sesh{4} = 1; Sesh{5} = [1,2]; Sesh{6} = [1 2];
                
            elseif strcmp(loc,'ac')
                TITLE = 'SOM + ArchT in AC';
                MNum = [3139 3144 3160 3163 3161 3162 3196];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [1]; Sesh{2} = [1]; Sesh{3} = [1,2]; Sesh{4} = 1; Sesh{5} = [1,2]; Sesh{6} = 2; Sesh{7} = [3,4];
                
            end
        elseif strcmp(mouseline,'camk2')
            if strcmp(loc,'ic')                       
                TITLE = 'CaMKII + ArchT in IC';
                MNum = [3153 3180 3181 3227 3228 3229];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 2; Sesh{2} = [1,2]; Sesh{3} = [1,2]; Sesh{4} = [1 2 3 4]; Sesh{5} = [1 2 3 4]; Sesh{6} = [1 2 3 4];
            elseif strcmp(loc,'ac')                
                TITLE = 'CaMKII + ArchT in AC';
                MNum = [3180 3181 3227 3228 3229];
                Sesh = cell(1,length(MNum));
                Sesh{1} = [3,4]; Sesh{2} = 3; Sesh{3} = [5]; Sesh{4} = [5 6]; Sesh{5} = [5 6];
            end
        elseif strcmp(mouseline,'pv')
            if strcmp(loc,'ic')                
                TITLE = 'PV + ArchT in IC';
                MNum = [3172 3182 3183 3209];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 1; Sesh{2} = [1 2]; Sesh{3} = [1 2]; Sesh{4} = [1 2];
            elseif strcmp(loc,'ac')
                TITLE = 'PV + ArchT in AC';
                MNum = [3172 3182 3183 3209];
                Sesh = cell(1,length(MNum));
                Sesh{1} = 2; Sesh{2} = [3,4]; Sesh{3} = 3; Sesh{4} = [3 4];

            end
        end
    end
    
end

