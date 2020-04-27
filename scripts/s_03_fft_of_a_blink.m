% inizializzo alcune variabili
clear 

close all

srate = 500;
 
% list so5me frequencies
frex = [ 5 ];

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [ 70 ];


jitter_factor= 3

phases = 2*pi*jitter_factor;

% define time...
time=-0.2:1/srate:1;

 % per simulare del segnale transiente decido la in cui sar? effettivamente
 % presente del segnale. Specifico inizio e fine della finestra (solo due
 % numeri).
 Effectwin=[0.4 0.9]
 
 % aggiungo noise alle waveform create
 noise_level=0
 
 % per plottare (o meno) i singoli trials
 plot_single_trials=0;

 
% custom ylim 
 myylim=[-10 10]


% Creo una matrice con tre sine waves, con tutti i parametri uguali, tranne
% la fase che ? definita dal vettore phases (vedi sopra).


sine_waves = amplit * sin(2*pi*frex.*time + phases);

timepoints_edges=dsearchn( time', Effectwin');
timepoints=timepoints_edges(1):timepoints_edges(2);
 
% Creo delle sine waves in cui setto a zero alcune parti, replicando le
% sine waves originali
sine_waves_eff=sine_waves;

% simula dei blink (mettendo 0 invece puoi vedere l'effetto di fft su
% segnale.
if 1
sine_waves_eff(: , setdiff(1:length(sine_waves_eff), timepoints))=0;

% per evitare discontinuit? metto anche a 0 eventuali valori negativi
sine_waves_eff(sine_waves_eff<0) = 0;
end

%aggiungo noise
% commenta le due righe seguenti per togliere il noise
sine_waves_eff=sine_waves_eff+rand(size(sine_waves_eff)).*noise_level; 
 

figure 
plot(time, sine_waves_eff)
set(gca, 'ylim', [-100, 100])
%title('Error','fontweight','normal');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 3])
print(['../Figures/s_03_blink'], '-djpeg', '-r300');

figure
% compute FFT
sig_fft=fft(sine_waves_eff);
hz = linspace(0,srate/2,floor(length(time)/2)+1);
plot(hz,abs(sig_fft(1:length(hz))*2))
set(gca, 'xlim', [0, 60])


%title('Error','fontweight','normal');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 3])
print(['../Figures/s_03_FFT_blink'], '-djpeg', '-r300');










