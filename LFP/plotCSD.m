function plotCSD(CSD,sigma,xrange,yrange)
%Plot current source density, both as line plots and heatmap
%
%Inputs:
%   CSD = each row is a tetrode/electrode arranged from top to bottom
%   sigma = sigma for imgausfilt smoothing. Value of 0 is no smoothing.
%   xrange = range on the xaxis (ms)
%   yrange = tetrode/electrode numbers
%

subplot(1,2,1);
offset = NaN(1,size(CSD,1));
offset(1) = 0;
for i = 1:size(CSD,1)
    hold on;
    tetID = size(CSD,1) + 1 - i; %Reverse the order so plots bottom of probe at bottom of plot
    if i > 1
        offset(i) = max(abs(CSD(tetID,:))) + max(abs(CSD(tetID+1,:))) + offset(i-1);
    end
    plot(xrange,CSD(size(CSD,1) + 1 - i,:) + offset(i))
end
line([0 0],[min(CSD(end,:)) max(CSD(1,:))+offset(end)],'Color','k','LineStyle','--')
set(gca,'TickDir','out')
xlabel('Time from tone onset (ms)')
axis tight

subplot(1,2,2);
if sigma > 0
    imagesc(xrange,yrange,imgaussfilt(CSD,sigma));
elseif sigma == 0
    imagesc(xrange,yrange,CSD)
end
line([0 0],[yrange(1)-0.5 yrange(end)+0.5],'Color','k','LineStyle','--')
set(gca,'TickDir','out')
MM = max(CSD(:));
mm = min(CSD(:));
lim = max(abs(MM),abs(mm));
caxis([-lim,lim])
ylabel('Tetrode')
xlabel('Time from tone onset (ms)')
box off
colorbar

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[1600 1000]);
set(gcf,'Position',[0 0 1600 1000]);