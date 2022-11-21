%% FFT script
clear all; close all; clc;

SampleTime = 1/2048;

load('data/20181023_1025.mat')
FRF_y1 = ans';
% 1. Time
% 2. enc0 [y1]
% 3. enc1 [y2]
% 4. control output (dac input [u]
% 5. error (control input) [e]
% 6. disturbance input (inserted at control input) [d]
% 7. input/reference [r]

%% TFestimate y1
nfft = 10*(1/SampleTime); % Window length

[S_estimate_y2,f] = tfestimate(FRF_y1(:,6), FRF_y1(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[S_coh_y2, ~] = mscohere(FRF_y1(:,6), FRF_y1(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

[PS_estimate_y2,~] = tfestimate(FRF_y1(:,6), FRF_y1(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_coh_y2, ~] = mscohere(FRF_y1(:,6), FRF_y1(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

H_estimate_y1 = (-PS_estimate_y2)./S_estimate_y2;
H_coh_y1 = S_coh_y2.*PS_coh_y2;

if (true)
    figure()
    subplot(3,1,1)
    semilogx(f,db(S_estimate_y2));
    title('Sensitivity y2')
    ylabel('Amplitude [db]')
    grid on;
    subplot(3,1,2)
    semilogx(f,(angle(S_estimate_y2)*180/pi));
    ylabel('Angle [deg]')
    grid on;
    subplot(3,1,3)
    semilogx(f,S_coh_y2);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(-PS_estimate_y2));
    grid on;
    title('Process sensitivity y2')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(-PS_estimate_y2)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,PS_coh_y2);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(H_estimate_y1));
    grid on;
    title('Plant y2')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(H_estimate_y1)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,H_coh_y1);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
end