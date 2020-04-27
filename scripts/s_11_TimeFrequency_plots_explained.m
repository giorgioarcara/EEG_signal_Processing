%%  This code is an adaptation of code by Mike Cohen  (chapter13.m, figure 13.11).
% all credits to him.
% all blames for mistakes to me.

close all
clear all

EEGdata = load('../example_data/data_S131_trial001.mat');
ChannelMat = load('../example_data/channel.mat');

% calculate sampling rate from points in TIme VEctor
srate = ((length(EEGdata.Time)-1).*1000)/ ( (-EEGdata.Time(1).*1000) + (EEGdata.Time(end).*1000 ));


close all

% definitions, selections...
chan2use = 'O1';

min_freq =  2;
max_freq = 80;
num_frex = 30;

% define wavelet parameters
time = -1:1/srate:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters

n_wavelet            = length(time);
n_data               = length(EEGdata.Time);
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
EEG_pnts = size(EEGdata.F,2);

eegfft = fft(reshape(EEGdata.F(strcmpi(chan2use,{ChannelMat.Channel.Name}),:,:),1,EEG_pnts),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG_pnts); % frequencies X time X trials

baseidx = dsearchn(EEGdata.Time',[-1.2 -1]');

% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG_pnts,1)).^2,2);
    %eegpower(fi,:) =
    %10*log10(temppower./mean(temppower(baseidx(1):baseidx(2)))); % decibel
    %change
    baseline = mean(temppower(baseidx(1):baseidx(2)));
    % eegpower(fi,:)=10*log10(temppower./baseline); % decibel change
    eegpower(fi,:) = ((temppower - baseline)/baseline.*100) ; %ERS/ERD
    %eegpower(fi,:) = temppower; %raw power
end;



% CREO GLI ASSI
% per usare la funzione surf, devo creare separatamente degli assi di x
% (tempo) e y (frequenze), in modo che siano fatti da numeri interi.
% per tale scopo uso dsearchn, che restituisce l'indice del valore pi? vicino 
% (nel vettore) al numero (ms o frequenza intera) che io specifico.

% creo asse delle x (tempo)
Timeindices=dsearchn(EEGdata.Time', (-1.5:0.5:1.5)');
Timelabels=(-1.5:0.5:1.5);

% creo asse delle y (frequenze)
frexindices=dsearchn(frex', [ 2 5 8 10 15 20 40 50 80]'); 
frexlabels=[ 2 5 8 10 15 20 40 50 80];


%% Raw EEG data (1 epoch)
figure
plot(EEGdata.Time, EEGdata.F(strcmpi(chan2use,{ChannelMat.Channel.Name}),:).*100000)
set(gca, 'XTick', Timelabels, 'YLim', [-5 5])
title('Single Trial')

pow_lims = 'auto';
%pow_lims = [-12000 12000];


if strcmp(pow_lims, 'auto')
    pow_max = max(abs(eegpower(:)));
    pow_lims = [-pow_max pow_max];
end;



%% FIGURE 2
figure
sel_freq = 10; % to show this frequency
sel_freq_ind = dsearchn(frex', sel_freq'); 

plot(EEGdata.Time, eegpower(sel_freq_ind,:))
set(gca, 'Xtick', Timelabels)
title(['Power at' num2str(sel_freq), 'hz'])


%% 
% GRAFICO DI UNA SOLA BANDA DI FREQUENZA

%% FIGURE 3
% creo duplicato eegpower in cui setto NaN
% tutte le freuqenze tranne una piccola banda
figure
eegpower_masked=eegpower;
sel_freq_stripe_ind = (sel_freq_ind-1):1:(sel_freq_ind+1);
eegpower_masked(setdiff(1:size(eegpower_masked, 1), sel_freq_stripe_ind),:)=NaN;
surf(EEGdata.Time, 1:size(eegpower_masked,1),eegpower_masked)
set(gca,'clim', pow_lims, 'zlim', pow_lims, 'YLim', [1, size(eegpower,1)],  'ytick', frexindices, 'yticklabel', frexlabels);
shading interp
colormap jet

%% FIGURE 4
% GRAFICO DI TUTTE LE FREQUENZE
figure
surf(EEGdata.Time, 1:size(eegpower,1),eegpower)
set(gca,'clim', pow_lims, 'zlim', pow_lims, 'YLim', [1, size(eegpower,1)], 'Ytick', frexindices, 'YTicklabel', frexlabels);
shading interp
colormap jet
rotate3d
set(gca,'fontsize',15)
view([-30, 40]);

%% FIGURE 5
% GRAFICO DI TUTTE LE FREQUENZE
figure
contourf(EEGdata.Time, 1:size(eegpower,1),eegpower)
set(gca,'clim', pow_lims, 'zlim', pow_lims, 'YLim', [1, size(eegpower,1)], 'Ytick', frexindices, 'YTicklabel', frexlabels);
shading interp
colormap jet
set(gca,'fontsize',15)
colorbar


%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 10])
%print('TF_3d', '-djpeg', '-r0');

