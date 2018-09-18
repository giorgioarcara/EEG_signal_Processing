% inizializzo alcune variabili
clear

srate = 100;

% list so5me frequencies
frex = 40;

% list some random amplitudes... make sure there are
% the same number of amplitudes as there are frequencies!
amplit = [ 1 ];

% phases... trial 1 e 2 sono in fase, il 3 ? sfasato di pi/8.
phases = [pi/8];

% define time...
time=-1:1/srate:1;



% Creo una matrice con tre sine waves, con tutti i parametri uguali, tranne
% la fase che ? definita dal vettore phases (vedi sopra).

sine_waves = zeros(length(frex),length(time)); % remember: always initialize!
for phi=1:length(phases)
    sine_waves (phi,:) = amplit * sin(2*pi*frex.*time + phases(phi));
end

% per simulare del segnale transiente decido la in cui sar? effettivamente
% presente del segnale. Specifico inizio e fine della finestra (solo due
% numeri).
Effectwin=[-1 1]


timepoints_edges=dsearchn( time', Effectwin');
timepoints=timepoints_edges(1):timepoints_edges(2);

% Creo delle sine waves in cui setto a zero alcune parti, replicando le
% sine waves originali
sine_waves_eff=sine_waves;
sine_waves_eff(: , setdiff(1:length(sine_waves_eff), timepoints))=0;

% aggiungo noise
% commenta le due righe seguenti per togliere il noise
noise_level=0
sine_waves_eff=sine_waves_eff+rand(size(sine_waves_eff)).*noise_level;

figure
plot(time, sine_waves_eff(:,:))
set(gca, 'Ylim', [-(amplit+2), +amplit+2])



% PROVA WAVELET
% parameters...
srate = srate; % sampling rate in Hz
f     = 40; % frequency of wavelet in Hz
time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate
s     = 4/(2*pi*f);

% and together they make a wavelet
wavelet = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2));

plot(real(wavelet))


kernel = wavelet;

TFresults=zeros(size(sine_waves_eff));
TFresults_conv = zeros(size(sine_waves_eff));


figure

%% LOOP DISABILITATION
%for i =1:size(TFresults,1)

i=1
% convolution with conv function
TFresults(i,:)=conv(sine_waves_eff(i,:), wavelet, 'same');



% data that we'll use for convolution (must be zero-padded).
dat4conv = [zeros(1,length(kernel)-1) sine_waves_eff(i,:) zeros(1,length(kernel)-1) ];



% used for cutting the result of convolution
half_of_kernel_size = ceil((length(kernel)-1)/2);

% initialize convolution output
convolution_result = zeros(1,length(sine_waves_eff(i,:))+length(kernel)-1);

%% plot data to convolve

subplot(3,1,1)
plot(1:length(dat4conv), dat4conv)
%%highlight_padding
padding = repmat(NaN, 1, length(dat4conv));
padding(1:length(kernel)-1)=0;
padding(end-length(kernel)-1:end)=0;
plot(1:length(dat4conv), dat4conv)
hold on
plot(1:length(dat4conv), padding, '.', 'col', 'black', 'MarkerSize', 6);
set(gca, 'xLim', [1, length(dat4conv)]);
hold off


% run convolution (note that kernel is flipped backwards)
for ti=1:length(convolution_result)
    
    kern_real = real(kernel);
    subplot(3,1,2)
    % plot kernel
    
    plot(ti:length(kernel)+ti-1, real(kernel(end:-1:1)), 'col', 'red');
    set(gca, 'xLim', [1, length(dat4conv)]);
    convolution_result(ti) = sum(dat4conv(ti:ti+length(kernel)-1).*kernel(end:-1:1));
    
    % plot convolution
    subplot(3,1,3);
    plot(1+half_of_kernel_size:ti+half_of_kernel_size, abs(convolution_result(1:ti)), '.', 'col', 'red', 'MarkerSize', 2);
    set(gca, 'xLim', [1, length(dat4conv)], 'YLim', [-10,100]);
    pause(0.005)
    
    
    
    % end
    
    TFresults_conv(i, :)=convolution_result(half_of_kernel_size+1:end-half_of_kernel_size);
    
    
    
    
end;


%figure
%plot(time, real(wavelet))


% POWER dell segnale dopo convoluzione con wavelet a frequenza pari al
% segnale.
%figure
%plot(time, abs(TFresults(1,:)).^2)



%figure
%plot(time, abs(TFresults_conv(1,:)).^2)

% il power medio non risentir? del problema degli ERP.


% NOTA: questo script ? un po' circolare. Il modo in cui ho generato il
% segnale ? (a meno dell'effetto di una gaussian window) quello utilizzato
% alla base delle wavelet.