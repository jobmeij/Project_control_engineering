%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%   4CM00 - Control Engineering                 %
%  Practicum - plot error and analyse the error %
%                                               %
%   Author: Marcel van Wensveen                 %
%       &   Job Meijer                          %
%   Date:   24-10-2019                          %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all hidden, clc

%% Load data

load('24 okt\Controller_23hz_withFF.mat')
Data_withFF = Data;
load('24 okt\Controller_23hz_withoutFF.mat')
Data_withoutFF = Data;

%load('MarcelControl_FF_NoNotch_No5HzControlNotch.mat')

load('ScanningSetpointStruct.mat')

%% Plot time domain error
%Data_withFF = Data;
StartTime = 0; % seconds
StopTime = 3.6; % seconds
figure(1)
plot(Data_withFF.Time, Data_withFF.e,Data_withoutFF.Time, Data_withoutFF.e, Setpoint.Time, Setpoint.Position/40)
xlim([StartTime StopTime])
xlabel('Time [sec]')
ylabel('Amplitude [Rad]')
title('Error during the scanning movement with the 23Hz controller')
legend('Error with controller 2 + FF','Error with improved controller + FF','Scaled position setpoint for reference', 'Location','best')

%% Plot time domain of setpoint and output y2
StartTime = 0; % seconds
StopTime = 3.6; % seconds
figure(1)
plot(Data_withFF.Time, Data_withFF.y2,Data_withoutFF.Time, Data_withoutFF.y2, Setpoint.Time, Setpoint.Position)
xlim([StartTime StopTime])
xlabel('Time [sec]')
ylabel('Amplitude [Rad]')
title('Error during the scanning movement with the 23Hz controller')
lege

%% Plot cummulative PSD of the error
Ts = Data.Time(2)-Data.Time(1);
fftDuration = 1.6; % seconds
N_fft = round(fftDuration/Ts);

[PSD_withFF,f] = pwelch(Data_withFF.e, hann(N_fft), N_fft/2, N_fft,1/Ts, 'Power');

Forward_PS_withFF = abs(PSD_withFF(1));
for i = 2:length(PSD_withFF)
    Forward_PS_withFF(i) = Forward_PS_withFF(i-1)+abs(PSD_withFF(i));
end
Forward_PS_withFF = sqrt(Forward_PS_withFF);
Backward_PSD_withFF(length(PSD_withFF)) = abs(PSD_withFF(end));
for i = length(PSD_withFF)-1:-1:1
    Backward_PSD_withFF(i) = Backward_PSD_withFF(i+1)+abs(PSD_withFF(i));
end
Backward_PSD_withFF = sqrt(Backward_PSD_withFF);

[PSD_withoutFF,f] = pwelch(Data_withoutFF.e, hann(N_fft), N_fft/2, N_fft,1/Ts, 'Power');

Forward_PS_withoutFF = abs(PSD_withoutFF(1));
for i = 2:length(PSD_withoutFF)
    Forward_PS_withoutFF(i) = Forward_PS_withoutFF(i-1)+abs(PSD_withoutFF(i));
end
Forward_PS_withoutFF = sqrt(Forward_PS_withoutFF);
Backward_PSD_withoutFF(length(PSD_withoutFF)) = abs(PSD_withoutFF(end));
for i = length(PSD_withoutFF)-1:-1:1
    Backward_PSD_withoutFF(i) = Backward_PSD_withoutFF(i+1)+abs(PSD_withoutFF(i));
end
Backward_PSD_withoutFF = sqrt(Backward_PSD_withoutFF);

figure(2)
subplot(3,1,1)
semilogx(f,db(PSD_withFF),f,db(PSD_withoutFF))
title('PSD of the error with controller of 10Hz')
legend('Controller 2 + FF','Improved controller + FF')
ylabel('Amplitude')
grid on;
subplot(3,1,2)
semilogx(f,Forward_PS_withFF,f,Forward_PS_withoutFF)
title('Forward commulative PSD')
ylabel('Amplitude')
grid on;
subplot(3,1,3)
semilogx(f,Backward_PSD_withFF,f,Backward_PSD_withoutFF)
title('Backward commulative PSD')
ylabel('Amplitude')
xlabel('Frequency [Hz]')
grid on;
