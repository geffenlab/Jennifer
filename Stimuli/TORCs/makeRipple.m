function [y,A,S,Y]=makeRipple(duration, fs,d,omega,w,ripple_phi,ntones,fmin,oct)
% This code creates moving ripples, which are periodic auditory gratings each
% having a spectro-temporal profile modulated sinusoidally in spectrum and in time (See Klein et al. 2000).
%
% Inputs:
%   duration = Length of ripples (in seconds)
%   fs = sample frequency (seconds)
%   d = modulation depth
%   ripple_phi = phase
%   w = angular frequency (Hz) [temporal modulation]
%   ntones = number of tones (essentially frequency resolution)
%   fmin = min frequency (Hz)
%   oct = number of octaves TORC spans

%% Set parameters
% duration = 0.25; % in seconds
% fs = 44100; % sample frequency
% 
% d =0.5; % modulation depth % of mean
% omega = 1.4; % cycles/octave spectral density
% ripple_phi = 0; % phase
% w = [-4,-8,-12,-16,-20]; % angular frequency in Hz
% 
% ntones = 500; % number of tones
% fmin = 250; % minimum frequency
% oct=7; % number of octaves TORC spans


fmax = fmin*2^oct; % max frequency based on number of octaves
t=linspace(0,duration,duration*fs); % time bins
f=exp(linspace(log(fmin),log(fmax),ntones)); % distribution of frequencies
gamma = rand(1, ntones);
% gamma=normrnd(0,1,[1,ntones]);
phi = 2*pi*rand(1,ntones); % Randomise phase
% phi = 2*pi*ones(1,ntones);
% phi= normrnd(0,2*pi,[1,ntones]);
T=repmat(t,[ntones,1]);
F=repmat(f,[length(t),1]);%.*T';
%F = np.tile(f, (t.size, 1)).T
Gamma=repmat(gamma,[length(t),1]);%.*T';
% Gamma = np.tile(gamma, (t.size, 1)).T
Phi=repmat(phi,[length(t),1]);%.*T';
% Phi = np.tile(phi, (t.size, 1)).T


X = log2(F / fmin);

% Create auditory gratings
A = 1 + d*sin(2 * pi * ((-w .* T') + omega.* X) + ripple_phi);
%     figure
%      imagesc(flipud(A'))

% Create the amplitude modulation of each frequency
S = Gamma.* sin(2 * pi * F.*T' + Phi) ./ sqrt(F);

% Apply amplitude modulation to gratings
Y = A'.*S';

y = sum(Y);
A=A';
% figure
% imagesc(Y)
% soundsc(y/10,fs)
% imagesc(A)