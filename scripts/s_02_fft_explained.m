clear 

close all

srate = 500;
 
% list some frequencies
frex = [1 10 40];

% list some random amplitudes... make sure there are 
% as many amplitude values as frequencies
amplit = [ 100 20 30 ];

phases = pi

% define time...
time=0:1/srate:1;


for iFreq=1:length(frex)
    sine_waves (iFreq,:) = amplit(iFreq) * sin(2*pi*frex(iFreq).*time + phases);
end

% plot composing waves
n_sine_waves = size(sine_waves,1);
summed_signal = sum(sine_waves, 1);


figure
for iSin = 1:((n_sine_waves))
    subplot(n_sine_waves+1, 1, iSin)
    plot(time, sine_waves(iSin,:));
    set(gca, 'Ylim',[-max(summed_signal), max(summed_signal)])
end;
subplot(n_sine_waves+1, 1, iSin+1)
plot(time, sum(sine_waves, 1), 'Color', 'red');
set(gca, 'Ylim',[-max(summed_signal), max(summed_signal)])
title('Sum')
print -dpng ../Figures/s_02_Figure1.png;



figure
sig_fft=fft(sum(sine_waves, 1));
hz = linspace(0,srate/2,floor(length(time)/2)+1);
plot(hz,abs(sig_fft(1:length(hz))*2))
set(gca, 'Xlim', [0, 60])
print -dpng ../Figures/s_02_Figure2.png;

