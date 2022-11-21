%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%   4CM00 - Control Engineering                 %
%   Practicum - Single controller simulation    %
%                                               %
%   Author: Marcel van Wensveen                 %
%       &   Job Meijer                          %
%   Date:   17-10-2019                          %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all hidden, clc

%% Create PSD of error

load('Data\MarcelControl_FF_WithNotch_No5HzControlNotch.mat')
Ts = Data.Time(2)-Data.Time(1);
Duration = 1.8; % seconds
N_fft = round(Duration/Ts);
%freqvec = linspace(1/Duration,1/(2*Ts),1/Ts)';
%error =Data.e(100000:end); 

[PSD,f] = pwelch(Data.e, hann(N_fft), N_fft/2, N_fft,1/Ts, 'Power');
%PSD = PSD;

Forward_PSD = abs(PSD(1));
for i = 2:length(PSD)
    Forward_PSD(i) = Forward_PSD(i-1)+abs(PSD(i));
end
Forward_PSD = sqrt(Forward_PSD);
Backward_PSD(length(PSD)) = abs(PSD(end));
for i = length(PSD)-1:-1:1
    Backward_PSD(i) = Backward_PSD(i+1)+abs(PSD(i));
end
Backward_PSD = sqrt(Backward_PSD)

% plot Cummulative PDS
figure(2)
subplot(3,1,1)
semilogx(f,db(PSD))
title('PSD of the error')
ylabel('Amplitude')
grid on;
subplot(3,1,2)
semilogx(f,Forward_PSD)
title('Forward commulative PSD')
ylabel('Amplitude')
grid on;
subplot(3,1,3)
semilogx(f,Backward_PSD)
title('Backward commulative PSD')
ylabel('Amplitude')
xlabel('Frequency [Hz]')
grid on;

%% plot error
figure(5)
plot(Data.Time, Data.e)


