%ICA explained
% inizializzo alcune variabili
close all
clear all

srate = 500;

n_sources = 2

% list a couple of frequencies (they will be used to create signals).
frex = [7 3];

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [ 10 30 ];

Ntrials= 50

jitter_factor= 6

phases = [ pi, pi/2 ]; % same as frex and amplit

% define time...
time=-1:1/srate:1;

 
% custom ylim 
 myylim=[-10 10]

 
% generate sine waves from parameters above.
 sine_waves = zeros(length(frex),length(time)); % remember: always initialize!
for n=1:n_sources
    sine_waves (n,:) = amplit(n) * sin(2*pi*frex(n).*time + phases(n));
end

my_ylim=[-50,50]

figure
subplot(121)
plot(time, sine_waves(1,:));
title('A', 'FontSize', 18)
set(gca, 'ylim', my_ylim);
subplot(122)
plot(time, sine_waves(2,:));
title('B', 'FontSize', 18)
set(gca, 'ylim', my_ylim);

% create mixed signal
S1 = 0.5 .* sine_waves(1,:) + 0.8 .* sine_waves(2,:)
S2 = 0.4 .* sine_waves(1,:) + 0.2 .* sine_waves(2,:)

S = [S1; S2];

figure
subplot(121)
plot(time, S1 )
title(' S1 = 0.5 * A + 0.8 * B', 'FontSize', 18)
set(gca, 'ylim', my_ylim);
subplot(122)
plot(time, S2 )
title(' S2 = 0.4 * A + 0.2 * B', 'FontSize', 18)
set(gca, 'ylim', my_ylim);


addpath('/Applications/fieldtrip-20161013/external/fastica/')

[icasig, A, W]  = fastica(S)

figure
subplot(121)
plot(time, icasig(1,:));
title('IC 1', 'FontSize', 18)
set(gca, 'ylim', my_ylim);
subplot(122)
plot(time, icasig(2,:));
title('IC 2', 'FontSize', 18)
set(gca, 'ylim', my_ylim);


% Creo una matrice con tre sine waves, con tutti i parametri uguali, tranne
% la fase che ? definita dal vettore phases (vedi sopra).


r = A*icasig;


figure
subplot(121)
plot(time, r(1,:));
title({'Sensor 1 (reconstructed)', 'A * ICs'}, 'FontSize', 14)
set(gca, 'ylim', my_ylim);
subplot(122)
plot(time, r(2,:));
title({'Sensor 2 (reconstructed)', 'A * ICs'}, 'FontSize', 14)
set(gca, 'ylim', my_ylim);


A_1 = A;
A_1(:,1) = 0;
r1 = A_1*icasig;


figure
subplot(121)
plot(time, r1(1,:));
title({'Sensor 1 without IC2', 'A1 * ICs'}, 'FontSize', 14)
set(gca, 'ylim', my_ylim);
subplot(122)
plot(time, r1(2,:));
title({'Sensor 2 without IC2', 'A1 * ICs'}, 'FontSize', 14)
set(gca, 'ylim', my_ylim);







