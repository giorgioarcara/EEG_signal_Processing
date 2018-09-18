% inizializzo alcune variabili
close all

srate = 500;
 
% list so5me frequencies
frex = [4];

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [10];

Ntrials= 10

jitter_factor= 10

phases = rand(1,Ntrials)*pi*jitter_factor;

% define time...
time=-0.2:1/srate:1;

 % per simulare del segnale transiente decido la in cui sar? effettivamente
 % presente del segnale. Specifico inizio e fine della finestra (solo due
 % numeri).
 Effectwin=[0.3 0.6]
 
 % aggiungo noise alle waveform create
 noise_level=1
 
 % per plottare (o meno) i singoli trials
 plot_single_trials=1;

 
% custom ylim 
 myylim=[-20 20]




% Creo una matrice con tre sine waves, con tutti i parametri uguali, tranne
% la fase che ? definita dal vettore phases (vedi sopra).

sine_waves = zeros(length(frex),length(time)); % remember: always initialize!
for phi=1:length(phases)
    sine_waves (phi,:) = amplit * sin(2*pi*frex.*time + phases(phi));
end


 
 
 timepoints_edges=dsearchn( time', Effectwin');
 timepoints=timepoints_edges(1):timepoints_edges(2);
 
% Creo delle sine waves in cui setto a zero alcune parti, replicando le
% sine waves originali
sine_waves_eff=sine_waves;
sine_waves_eff(: , setdiff(1:length(sine_waves_eff), timepoints))=0;

%aggiungo noise
% commenxta le due righe seguenti per togliere il noise
sine_waves_eff=sine_waves_eff+rand(size(sine_waves_eff)).*noise_level;





%% PLOT FIGURE


if plot_single_trials % plot single trials and ERP
    figure
    for i=1:Ntrials
        subplot(Ntrials+1, 1, i); % il +1 ? per poter plottare l'ERP sotto
        plot(time, sine_waves_eff(i,:), 'LineWidth', 2);
    end;

    % plot ERP
    subplot(Ntrials+1,1,i+1)

    plot(time, mean(sine_waves_eff, 1), 'red'); 
    title('ERP');
    set(gca, 'ylim', myylim);
else % plot just ERP
    figure
    ERP=mean(sine_waves_eff, 1);
    plot(time,ERP , 'red'); 
    title('ERP');
    set(gca, 'ylim', myylim);
    %set(gca, 'ylim', [min(ERP) max(ERP)]);
    set(gca, 'ylim', myylim);
end,



