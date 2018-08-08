function h = plotSTA(time,freq,STA,kwidth,lims)

% function h = plotSTA(x,y,STA,kwidth,lims)
% 
% Plots an STRF using x frequency bins, y time bins and intensity
% values in 2D matrix STA
% [optional]: Gaussian smoothing with kwidth kernel and limit on
% plot magnitudes between two numbers in lims


if ~exist('kwidth','var')
    kwidth = 0.1;
end
smoothSTA = imgaussfilt(STA,kwidth);

if ~exist('lims','var')
    lims = [min(smoothSTA(:)) max(smoothSTA(:))];
end

%kernel = fspecial('gaussian',kwidth,1);
%smoothSTA = imfilter(STA,kernel,'replicate');


h = surf(time,fliplr(freq),smoothSTA);
set(gca,'YScale','log');
set(h,'linestyle','none')

xlabel('Time (ms)');
ylabel('Frequency (kHz)');

axis tight;
view(2)

caxis(lims);
colormap(jet);