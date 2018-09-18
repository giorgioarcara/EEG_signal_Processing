%% figure 13.11
%cd C:\Users\Mamma\Desktop\Giorgio\AnalyzingNeuralTimeSeriesData_MatlabCode\code
cd('/Users/giorgioarcara/Documents/Da PC Fabio/AnalyzingNeuralTimeSeriesData_MatlabCode/code')
load sampleEEGdata.mat


% definitions, selections...
chan2use = 'o1';

min_freq =  2;
max_freq = 80;
num_frex = 30;

% define wavelet parameters
time = -1:1/EEG.srate:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials

baseidx = dsearchn(EEG.times',[-500 -200]');

% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end;



% CREO GLI ASSI
% per usare la funzione surf, devo creare separatamente degli assi di x
% (tempo) e y (frequenze), in modo che siano fatti da numeri interi.
% per tale scopo uso dsearchn, che restituisce l'indice del valore pi? vicino 
% (nel vettore) al numero (ms o frequenza intera) che io specifico.

% creo asse delle x (tempo)
Timeindices=dsearchn(EEG.times', (-1000:200:1500)');
Timelabels=(-1000:200:1500);

% creo asse delle y (frequenze)
frexindices=dsearchn(frex', [ 2 5 8 10 15 20 40 50 80]'); 
frexlabels=[ 2 5 8 10 15 20 40 50 80];

% GRAFICO DI UNA SOLA FREQUENZA 

% creo duplicato eegpower in cui setto NaN
% tutte le freuqenze tranne una piccola banda
plot(eegpower(8,:))
set(gca, 'xtick', Timeindices, 'xticklabel', Timelabels, 'ylim', [-6 6])
title('Power at 8 hz')





%% 
% GRAFICO DI UNA SOLA BANDA DI FREQUENZA

% creo duplicato eegpower in cui setto NaN
% tutte le freuqenze tranne una piccola banda
eegpower_masked=eegpower;
eegpower_masked([1:6 10:30],:)=NaN;
surf(1:size(eegpower_masked,2), 1:size(eegpower_masked,1),eegpower_masked)
set(gca,'clim', [-6 6], 'zlim', [-6 6], 'ylim', [1, size(eegpower,1)], 'xtick', Timeindices, 'xticklabel', Timelabels, 'ytick', frexindices, 'yticklabel', frexlabels);
shading interp
colormap jet


% GRAFICO DI TUTTE LE FREQUENZE
figure
surf(1:size(eegpower,2), 1:size(eegpower,1),eegpower)
set(gca,'clim', [-6 6], 'zlim', [-6 6], 'ylim', [1, size(eegpower,1)], 'xtick', Timeindices(1:2:end), 'xticklabel', Timelabels(1:2:end), 'ytick', frexindices([1,5,8]), 'yticklabel', frexlabels([1,5,8]));
shading interp
colormap jet
rotate3d
set(gca,'fontsize',15)
view([-30, 40]);


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 10])
print('TF_3d', '-djpeg', '-r0');

