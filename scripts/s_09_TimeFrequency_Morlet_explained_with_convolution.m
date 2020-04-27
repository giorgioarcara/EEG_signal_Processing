% inizializzo alcune variabili
% this script is largely inspired to the book "Analyzing Neural Time Series
% Data " by Mike Cohen.
clear

close all

srate = 100;

% list so5me frequencies
frex = 10;

% list some random amplitudes... make sure there are
% the same number of amplitudes as there are frequencies!
amplit = [ 0.5 ];

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
Effectwin=[-0.5 0.5]
timepoints_edges=dsearchn( time', Effectwin');
timepoints=timepoints_edges(1):timepoints_edges(2);

Effectsize=10

% Creo delle sine waves in cui setto a zero alcune parti, replicando le
% sine waves originali
sine_waves_eff=sine_waves;
sine_waves_eff(: , intersect(1:length(sine_waves_eff), timepoints))=Effectsize.*sine_waves_eff(: , intersect(1:length(sine_waves_eff), timepoints));


%  noise

noise_level=0
sine_waves_eff_n=sine_waves_eff+rand(size(sine_waves_eff)).*noise_level;

figure
subplot(2,1,1)
plot(time, sine_waves_eff(:,:))
set(gca, 'Ylim', [-max(abs(sine_waves_eff_n)), max(abs(sine_waves_eff_n))])
title('clean')
subplot(2,1,2)
plot(time, sine_waves_eff_n(:,:))
set(gca, 'Ylim', [-(max(abs(sine_waves_eff_n))), max(abs(sine_waves_eff_n))])
title('plus noise')




%  WAVELET
% parameters...
srate = srate; % sampling rate in Hz
f     = 10; % frequency of wavelet in Hz
time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate
s     = 6/(2*pi*f);

% and together they make a wavelet
wavelet = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2));

% To plot the wavelet
% plot(real(wavelet))


kernel = wavelet;

TFresults=zeros(size(sine_waves_eff_n));
TFresults_conv = zeros(size(sine_waves_eff_n));


my_t_step = 1; % set step for plotting figure
autoplay = 1; % set 1 if you want the plot goes automatically (no need of keypress).
create_gif = 0;
create_png = 0;

%% LOOP DISABILITATION
for i = 1:size(TFresults,1)
    
    % convolution with conv function (the code below is a repetition)
    TFresults(i,:)=conv(sine_waves_eff_n(i,:), wavelet, 'same');
   
    % data that we'll use for convolution (must be zero-padded).
    dat4conv = [zeros(1,length(kernel)-1) sine_waves_eff_n(i,:) zeros(1,length(kernel)-1) ];
    
    % used for cutting the result of convolution
    half_of_kernel_size = ceil((length(kernel)-1)/2);
    
    % initialize convolution output
    % NOTE: normally I would initialize with zeros. Here I initialize with
    % NaN for plotting reasons (when I plot the convolved results I want to
    % plot only results that were actually convolved, and some initiaized
    % zeros could be present, and plotted, if my_t_step is set not to 1.)
    
    convolution_result = NaN(1,length(sine_waves_eff_n(i,:))+length(kernel)-1);
    
    %% plot data to convolve
    h = figure
    n = 1 % I use this for gif
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
    set(gcf, 'color', 'white');
    hold off
    
    
    % run convolution (note that kernel is flipped backwards)
    for ti=1:my_t_step:length(convolution_result)
        
        kern_real = real(kernel);
        subplot(3,1,2)
        % plot kernel
        
        plot(ti:length(kernel)+ti-1, real(kernel(end:-1:1)), 'col', 'red');
        set(gca, 'xLim', [1, length(dat4conv)]);
        convolution_result(ti) = sum(dat4conv(ti:ti+length(kernel)-1).*kernel(end:-1:1));
        
        % plot convolution
        subplot(3,1,3);
        plot(1+half_of_kernel_size:ti+half_of_kernel_size, abs(convolution_result(1:ti)), '.', 'col', 'red', 'MarkerSize', 6);
        set(gca, 'xLim', [1, length(dat4conv)], 'YLim', [-10,250]);
        
        if autoplay
            pause(0.005);
        else
            pause
        end
        
        
        % end
        
        TFresults_conv(i, :)=convolution_result(half_of_kernel_size+1:end-half_of_kernel_size);
        
        if create_gif
            
            filename='../Figures/Wavelet_convolution.gif';
            
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if n == 1
                imwrite(imind,cm,filename,'gif', 'DelayTime',0.1, 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'WriteMode','append');
            end
            
        end;
        if create_png
            if n==1
                out_dir = '../Figures/Wavelet_conv_files/';
                mkdir(out_dir);
            end;
            filename=[out_dir, 'Wavelet_convolution', num2str(n), '.png'];
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'png');
        end
        
        n = n+1
        
    end;
    
    
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
if 0
    
    % controlla la relazione tra parte reale e immaginaria delle wavelet
    plot3(time, real(wavelet), imag(wavelet))
    
    
    plot3(time, real(wavelet), real(wavelet))
    
    % parte reale e immaginaria di wavelet sono massimamente dissimili
    sum(real(wavelet).*imag(wavelet))
    
    %diversamente reale e immaginaria sono simili a se stessi
    sum(imag(wavelet).*imag(wavelet))
    sum(real(wavelet).*real(wavelet))
    
end



