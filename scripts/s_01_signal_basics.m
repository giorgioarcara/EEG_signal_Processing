% inizializzo alcune variabili
clear 

close all

srate = 500;
 
% list some frequencies
frex = [10 10];

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [ 100 50 ];

phases = [0 pi]

% define time...
time=0:1/srate:1;


for iFreq=1:length(frex)
    sine_waves (iFreq,:) = amplit(iFreq) * sin(2*pi*frex(iFreq).*time + phases(iFreq));
end


myylim=[-100 100];

%% plot in separate subplots
figure
for iPlot=1:size(sine_waves, 1);
    subplot(size(sine_waves,1), 1, iPlot);
    plot(time, sine_waves(iPlot,:));
   set(gca, 'ylim', myylim,'fontsize',15);

end;


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 10])
print('../Figures/s_01_sine_waves_subplot', '-djpeg', '-r100');


%% overaly
figure
hold on
for iPlot=1:size(sine_waves, 1);
    plot(time, sine_waves(iPlot,:), 'LineWidth', 3);
end;
hold off
set(gca, 'ylim', myylim,'fontsize',15);


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
print('../Figures/s_01_sine_waves_overlay', '-djpeg', '-r100');


