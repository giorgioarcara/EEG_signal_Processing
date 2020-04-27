% genero un segnale che non ? sinusoidale, ma a "sawtooth".

T = 40*(1/40);
Fs = 1000;
dt = 1/Fs;
t = 0:dt:T-dt;
x = sawtooth(2*pi*50*t);


%Nota che lo spettro mostra molte armoniche. Le armoniche sono conseguenza 
% del fatto che il segnale non ? sinusoidale. Non so se la corrente di rete
% ? a sawtooth (o semplicemente contiene anche altre componenti che lo
% rendono non sinusoidale).

% L'aspetto chiave ? che un segnale non sinusoidale genera armoniche ad una
% FFT.

% Trovato qui https://en.wikipedia.org/wiki/Spectral_density

figure
subplot(1,2,1)
plot(t, x)

x_fft=abs(fft(x));
subplot(1,2,2)
plot(x_fft(1:length(x_fft)/2))