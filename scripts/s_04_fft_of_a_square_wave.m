% inizializzo alcune variabili
clear

close all

srate = 500;

% time vector
t = linspace(-pi,2*pi,121);
t_lin = t/pi;

% ampl vector
square_wave = 1.15*square(2*t);
sine_wave = 1.15*sin(2*t)


%% figure
figure
subplot(221)
%plot(t_lin, square_wave,'.-')
%hold
plot(t_lin, sine_wave) % sin wave
xlabel('t')
grid on

subplot(222)
% compute FFT
sig_fft=fft(sine_wave);
hz = linspace(0,srate/2,floor(length(t_lin)/2)+1);
plot(hz,abs(sig_fft(1:length(hz))*2), 'red')


subplot(223)
plot(t_lin, square_wave,'.-')
%hold
%plot(t_lin, sine_wave) % sin wave
xlabel('t')
%xlabel('t / \pi')
grid on

subplot(224)
% compute FFT
sig_fft=fft(square_wave);
hz = linspace(0,srate/2,floor(length(t_lin)/2)+1);
plot(hz,abs(sig_fft(1:length(hz))*2))


%title('Error','fontweight','normal');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print(['../Figures/s_04_square_wave_FFT'], '-djpeg', '-r300');










