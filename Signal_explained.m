% inizializzo alcune variabili
clear 

close all

srate = 500;
 
% list so5me frequencies
frex = [1 10 40];

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [ 100 100 30 ];

Ntrials= 50

jitter_factor= 6

phases = rand(1,Ntrials)*pi*jitter_factor;
phi = 1

% define time...
time=0:1/srate:1;

 % per simulare del segnale transiente decido la in cui sar? effettivamente
 % presente del segnale. Specifico inizio e fine della finestra (solo due
 % numeri).
 Effectwin=[0.3 0.5]
 
 % aggiungo noise alle waveform create
 noise_level=1
 
 % per plottare (o meno) i singoli trials
 plot_single_trials=0;

 
% custom ylim 
 myylim=[-100 100]


for iFreq=1:length(frex)
    sine_waves (iFreq,:) = amplit(iFreq) * sin(2*pi*frex(iFreq).*time + phases(phi));
end

plot(time, sum(sine_waves, 1));
set(gca, 'ylim', myylim)


% compute FFT
sig_fft=fft(sum(sine_waves, 1));
hz = linspace(0,srate/2,floor(length(time)/2)+1);
plot(hz,abs(sig_fft(1:length(hz))*2))


set(gca, 'xlim', [0, 60])







