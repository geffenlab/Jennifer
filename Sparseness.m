function S = Sparseness(mF,nStim);
%Sparseness measure from Weliky et al. (2003) [Take value between 0 and 1]. Higher values (closer to 1)
%indicate neuron responds to a narrower range of stimuli.

% INPUTS:
%   mF = vector of mean firing rate of a single neuron to each of the nStim
%        stimuli.
%   nStim = number of stimuli

% OUTPUTS:
%   S = Sparseness index


S = 1 - ((sum((mF./nStim)).^2)./(sum((mF.^2)./nStim)));