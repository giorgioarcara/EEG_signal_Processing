%% generate some signal

% set sampling rate
srate = 500;
 
% calculate Nyquist Frequency
Nyquist= srate/2

% list a frequency
frex = 80;

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [ 3 ];

% phases... trial 1 e 2 sono in fase, il 3 ? sfasato di pi/8.
phases = [pi/8]

% define time...
time=-1:1/srate:1;



% generate a sine_wave

sine_wave = amplit * sin(2*pi*frex.*time + phases);


%compute FFT
sine_fft=fft(sine_wave);

% determine x-axis of fft plot.
frequencies = linspace(0,Nyquist,length(sine_wave)/2+1);


figure
plot(frequencies, abs(sine_fft(1:length(frequencies))));
xlabel('Frequency (Hz)')


% resample data by selecting a subset of points
resample_factor=10

time_res=time(1:resample_factor:end);


sine_wave_res=sine_wave(1: resample_factor :end);

srate_res=srate/resample_factor

Nyquist_res=srate_res/2

sine_res_fft=fft(sine_wave_res);

frequencies_res = linspace(0, Nyquist_res, length(sine_wave_res)/2+1);
figure
plot(frequencies_res, abs(sine_res_fft(1:length(frequencies_res))));
xlabel('Frequency (Hz)')
title('Resampled data')


time2plot=-1:1/srate:-0.8
points2plot_or=dsearchn(time', time2plot');
points2plot_res=dsearchn(time_res', time2plot');

figure
subplot(3,1,1)
hold on
plot(time(points2plot_or), sine_wave(points2plot_or));
plot(time_res(points2plot_res), sine_wave_res(points2plot_res), 'o-', 'color', 'red', 'LineWidth', 2)
title('blue = original signal; red = resampled signal')
hold off

subplot(3,1,2)
plot(frequencies, abs(sine_fft(1:length(frequencies))));
xlabel('Frequency (Hz)')


subplot(3,1,3)
% NOTE: to make a comparison with the original data (and to show the
% aliasing) I set the same xlim of the original data.
plot(frequencies_res, abs(sine_res_fft(1:length(frequencies_res))));
set(gca, 'xlim', [0 Nyquist+1])
xlabel('Frequency (Hz) - Resampled data')








