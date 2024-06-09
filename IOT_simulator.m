[waveform,Fs] = generate_NBiot_UL();
%t = linspace(0,0.5*10^-3 , Fs);
%t = 0:1/Fs:0.08-1/Fs;
%waveform = upsample(waveform,3);
%Fs= Fs*4;
%t =  0:1/Fs:0.08-1/Fs;
%waveform = waveform' .* cos(2*pi*fc.*t);
%another approach for filtering is to use these functions 
%h = fir1(50 , [0.01 0.3]);
%freqz(h,100);
%waveform = filter(h,1,waveform);
%waveform = lowpass(waveform,0.33);

%transmitter starts here
upsampling_facotr = 5;
Fs= Fs*upsampling_facotr;
t =  0:1/Fs:0.08-1/Fs;
baseFc = 1920*10^3;
bandwidth = 200000;
%Fs = Fs * 4;
waveforms = [];
SNRs = [];
for i = 0:5
    [waveform,Fs] = generate_NBiot_UL();
    SNRs = [SNRs abs(mean(waveform.^2))];
    waveform = upsample(waveform,upsampling_facotr);
    Fs= Fs*upsampling_facotr;
    waveform = lowpass(waveform,0.2);
    waveform = waveform' .* cos(2*pi*(baseFc + i*bandwidth).*t);
    waveforms = [waveforms; waveform];
end

for i = 1:6
[signal,f] = pwelch(waveforms(i,:),2048,1024,2048,Fs);
signal = signal/max(abs(signal));
%figure;
plot(f/1e3,20*log10(signal))
hold on;
end

hold off;


signal_transmitted = sum(waveforms);
figure
[signal,f] = pwelch(signal_transmitted,2048,1024,2048,Fs);
signal = signal/max(abs(signal));
plot(f/1e3,20*log10(signal))

%%
%base station

% %testing here which filter is the best
% Hd = load("coffcs_LS200.mat","Hd");
% Hd = Hd.Hd;
% signal_recived = signal_transmitted .* cos(2*pi*(baseFc).*t);
% corr_factor = max(signal_transmitted)/max(signal_recived)
% signal_recived_filter = filter(Hd,1,signal_recived); %* 7.5373;
% %9.4116 -17.2091i
% %Noise_sig_filter = signal_recived - signal_transmitted;
% 
% snr_ls_200 = snr(signal_transmitted,ones(1,length(signal_recived_filter)))
% snr_org = snr(signal_recived,ones(1,length(signal_recived_filter)))
% %plot(abs(Noise_sig_filter))
% 
% 
% % [signal,f] = pwelch(signal_recived_filter,2048,1024,2048,Fs);
% % signal = signal/max(abs(signal));
% % plot(f/1e3,20*log10(signal))
% 
% [test_wave,~] = generate_NBiot_UL();
% 
% signal_downSampled = downsample(signal_recived_filter,4);
% 
% % figure
% % plot(abs(signal_downSampled))
% % hold on
% % plot(abs(test_wave(25:end,1))*(1))
% % legend
% 
% %figure
% [signal,f] = pwelch(signal_downSampled,2048,1024,2048,Fs);
% signal = signal/max(abs(signal));
% plot(f/1e3,20*log10(signal))
% hold on
% [signal,f] = pwelch(test_wave,2048,1024,2048,Fs);
% signal = signal/max(abs(signal));
% plot(f/1e3,20*log10(signal))
% 
% 
% noise_power = abs(mean(signal_recived_filter.^2)) - SNRs(1);
% snr_t = -abs(10*log10(SNRs(1)/noise_power))
% 
% fft_org = fft(test_wave,2048);
% fft_new = fft(signal_downSampled,2048);
% 
% snr_another_att = -10 * log10(abs(sum(fft_org(1:200).^2))/abs(sum(fft_new(1:200).^2)))


%actual reciever

signals_final = [];
Hd = load("coffcs_LS200.mat","Hd");
Hd = Hd.Hd;
signal_recived = signal_transmitted .* cos(2*pi*(baseFc).*t);
figure
[signal,f] = pwelch(signal_recived,2048,1024,2048,Fs);
signal = signal/max(abs(signal));
plot(f/1e3,20*log10(signal))
%filtering each device
for i = 1:6
    signal_enter_filter = signal_recived .* cos(2*pi*i*bandwidth.*t);
    signal_recived_filter = filter(Hd,1,signal_recived);
    %plot(abs(signal_recived_filter));
    %tesging something in the next 3 lines
    [signal,f] = pwelch(signal_recived_filter,2048,1024,2048,Fs);
    signal = signal/max(abs(signal));
    plot(f/1e3,20*log10(signal))
    hold on

    signal_downSampled = downsample(signal_recived_filter,4);
    
    signals_final = [signals_final; signal_downSampled];
end