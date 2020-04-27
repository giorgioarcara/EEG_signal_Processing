% inizializzo alcune variabili
close all
clear all

srate = 500;

% list so5me frequencies
frex = [10];

% list some random amplitudes... make sure there are
% the same number of amplitudes as there are frequencies!
amplit = [1];

Ntrials= 5

% jitter factor in ms
jitter_factor= 100

% define time...
time=-1:1/srate:3;

time_jitter = normrnd( 0, jitter_factor, [1, Ntrials]) ./1000;


% per simulare del segnale transiente decido la in cui sar? effettivamente
% presente del segnale. Specifico inizio e fine della finestra (solo due
% numeri).
Effectwin=[1 1.5]

Effectwin_jitter_factor = 10; % in ms


% aggiungo noise alle waveform create
noise_level=2
% per plottare (o meno) i singoli trials
plot_single_trials=1;



phases = 2*pi;
% Creo una matrice con le sine waves, con tutti i parametri uguali, tranne
% la fase che ? definita dal vettore phases (vedi sopra).

sine_waves = zeros(length(frex),length(time)); % remember: always initialize!
for tj=1:length(time_jitter)
    sine_waves (tj,:) = amplit * sin(2*pi*frex.*(time+time_jitter(tj))+ phases);
end


% Creo delle sine waves in cui setto a zero alcune parti, replicando le
% sine waves originali
sine_waves_eff=sine_waves;

eff_win = 'gaussian'
%eff_win = 'rect'


%% create effect windows (with jitter)
for iS = 1:size(sine_waves_eff)
    
    curr_Effectwin = Effectwin + normrnd( 0, Effectwin_jitter_factor, 1) ./1000; ; % in ms
    curr_Effect_ind = dsearchn( time', curr_Effectwin');
    curr_Effect_vec = curr_Effect_ind(1):curr_Effect_ind(2)
    
    if strcmp(eff_win, 'gaussian')
        sine_waves_eff(iS, curr_Effect_vec) = gausswin(length(curr_Effect_vec))'.* sine_waves_eff(iS, curr_Effect_vec);
        
    elseif strcmp(eff_win, 'rect')
        sine_waves_eff(iS, curr_Effect_vec) = rectwin(length(curr_Effect_vec))'.* sine_waves_eff(iS, curr_Effect_vec);
    end;
    sine_waves_eff(iS, setdiff(1:size(sine_waves_eff,2), curr_Effect_vec)) = 0; % set all other values to 0
   
end;


%% aggiungo noise
sine_waves_eff =sine_waves_eff + noise_level .* rand(size(sine_waves_eff));

% perform baseline correction
baseline = [-0.2, 0]
baseline_tp = dsearchn( time', baseline');

for iS = 1:size(sine_waves_eff, 1)
    sine_waves_eff(iS,:) = sine_waves_eff(iS,:) - mean(sine_waves_eff(iS, baseline_tp(1):baseline_tp(2)), 2);
end;




%% PLOT FIGURE ERP
% custom ylim
if noise_level==0
    myylim=[-max(amplit) max(amplit)]
else
    myylim = [-max(abs(sine_waves_eff(:))+1) max(abs(sine_waves_eff(:))+1)];
end;



if plot_single_trials % plot single trials and ERP
    figure
    for i=1:Ntrials
        subplot(Ntrials+1, 1, i); % il +1 ? per poter plottare l'ERP sotto
        plot(time, sine_waves_eff(i,:), 'LineWidth', 2);
        set(gca, 'Ylim', myylim);
    end;
    
    % plot ERP
    subplot(Ntrials+1,1,i+1)
    
    plot(time, mean(sine_waves_eff, 1), 'red');
    title('ERP');
    set(gca, 'Ylim', myylim);
else % plot just ERP
    figure
    ERP=mean(sine_waves_eff, 1);
    plot(time,ERP , 'red');
    set(gca, 'Ylim', myylim);
    title('ERP');
    set(gca, 'Ylim', myylim);
end,


%% TF

plot_tf = 1

if plot_tf
    %  WAVELET
    % parameters...
    srate = srate; % sampling rate in Hz
    f     = 5; % frequency of wavelet in Hz
    w_time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate
    s     = 18/(2*pi*f);
    
    % and together they make a wavelet
    wavelet = exp(2*pi*1i*f.*w_time) .* exp(-w_time.^2./(2*s^2));
    
    % To plot the wavelet
    % plot(real(wavelet))
    
    
    kernel = wavelet;
    
    TFresults=zeros(size(sine_waves_eff));
    
    for iS = 1:size(sine_waves_eff,1)
        TFresults(iS,:)=conv(sine_waves_eff(iS,:), wavelet, 'same');
    end;
    
    if plot_single_trials % plot single trials and ERP
        figure
        for i=1:Ntrials
            subplot(Ntrials+1, 1, i); % il +1 ? per poter plottare l'ERP sotto
            plot(time, abs(TFresults(i,:)), 'LineWidth', 2);
            set(gca, 'Ylim', [0 100]);
        end;
        
        % plot Average
        subplot(Ntrials+1,1,i+1)
        
        plot(time, mean(abs(TFresults), 1), 'red');
        set(gca, 'Ylim', [0 100]);
        title('TF average');
        %set(gca, 'ylim', myylim);
    else % plot just ERP
        figure
        TF=mean(sine_waves_eff, 1);
        plot(time,TF , 'red');
        %set(gca, 'Ylim', myylim);
        title('TF average');
    end,
end;
