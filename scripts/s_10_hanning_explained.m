clear all
close all
load('../example_data/data_S131_trial001.mat')

data2plot =  F(28,:);

figure

xlim = [Time(1), Time(end)];

subplot(3,1,1)
plot(Time, data2plot)
set(gca,'fontsize',12, 'Color', 'white', 'Xlim', xlim)

subplot(3,1,2)
plot(Time, hanning(length(data2plot))')
set(gca,'fontsize',12, 'Color', 'white' ,'Xlim', xlim)


subplot(3,1,3)
plot(Time, hanning(length(data2plot))' .* data2plot)
set(gca,'fontsize',12, 'Color', 'white',  'Xlim', xlim)

set(gcf,'Color', 'white', 'PaperUnits','inches','PaperPosition',[0 0 6 6])
print('../Figures/hanning_explained', '-dpng', '-r190');


