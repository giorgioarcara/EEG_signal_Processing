% WAVELET
% parameters...
srate = 500; % sampling rate in Hz
f     = 5; % frequency of wavelet in Hz
time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate
s     = 6/(2*pi*f);

% and together they make a wavelet
signal = exp(2*pi*1i*f.*time);
gaussian_win = exp(-time.^2./(2*s^2));

wavelet = signal .* gaussian_win; 

figure
subplot(3,1,1)
plot(time, real(signal))
set(gca, 'FontSize', 15)

subplot(3,1,2)
plot(time, gaussian_win)
set(gca, 'FontSize', 15)

subplot(3,1,3)
plot(time, wavelet)
set(gca, 'FontSize', 15)


print('wavelet', '-djpeg', '-r100');


%% plot how wavelet are really (with real and imaginary part).
plot3(time, real(wavelet), imag(wavelet))
rotate3d

