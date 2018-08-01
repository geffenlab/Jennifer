function [MNum,Sesh,TITLE] = loadsessions2(opsin, loc)
%Load mouse and session numbers for different experiment conditions for
%retroAAV + opsin experiments
%
%Inputs: (ALL LOWERCASE STRINGS)
%   opsin: name of opsin, 'chr2' or 'archt'
%   loc: recording location, 'ac' or 'ic
%
%Outputs:
%   MNum: list of mouse IDs (1xN array)
%   Sesh: list of sessions for each mouse ID (1xN cell array)
%   TITLE: string listing the input experiment info
%

if strcmp(opsin,'chr2')
        if strcmp(loc,'ic')
            TITLE = 'Feedback ChR2 in IC';
            MNum = [3222 3223 3224 3225 3226];
            Sesh = cell(1,length(MNum));
            Sesh{1} = [1 2 3 4]; Sesh{2} = [1 2 3 4]; Sesh{3} = [1 2 3 4]; Sesh{4} = [1 2 3 4]; Sesh{5} = [1 2 3 4]; 
        elseif strcmp(loc,'ac')
            TITLE = 'Feedback ChR2 in AC';
            MNum = [3122 3126];
            Sesh = cell(1,length(MNum));
            Sesh{1} = [5]; Sesh{2} = [5];
        end

elseif strcmp(opsin,'archt') %ArchT
        if strcmp(loc,'ic')
            TITLE = 'Feedback ArchT in IC';
            MNum = [3230 3231 3232];
            Sesh = cell(1,length(MNum));
            Sesh{1} = [1 2 3 4]; Sesh{2} = [1 2 3 4]; Sesh{3} = [1 2 3 4];
        elseif strcmp(loc,'ac')
            TITLE = 'Feedback ArchT in AC';
            MNum = [3130 3131 3132];
            Sesh = cell(1,length(MNum));
            Sesh{1} = [5]; Sesh{2} = [5]; Sesh{3} = [5];
        end
 end