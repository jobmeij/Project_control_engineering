%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
%   4CM00 - Control Engineering      %
%   Practicum - Model estimation     %
%                                    %
%   Author: Marcel van Wensveen      %
%       &   Job Meijer               %
%   Date: 18-10-2019                 %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load measured data and estimate the FRF

clear all; close all; clc;

% Init data
SampleRate = 2048;
SampleTime = 1/SampleRate;

load('data/20181012_1538.mat')
FRF_y1_control = ans';
load('data/20181012_1544.mat')
FRF_y2_control = ans';
clear ans;

nfft = 16*(1/SampleTime); % Window length

% system 1 (motor to mass 1)
[S_estimate_y1c,f] = tfestimate(FRF_y1_control(:,6), FRF_y1_control(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[S_coh_y1c, ~] = mscohere(FRF_y1_control(:,6), FRF_y1_control(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_estimate_y1c,~] = tfestimate(FRF_y1_control(:,6), FRF_y1_control(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_coh_y1c, ~] = mscohere(FRF_y1_control(:,6), FRF_y1_control(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
H_estimate_y1c = (-PS_estimate_y1c)./S_estimate_y1c;
H_coh_y1c = S_coh_y1c.*PS_coh_y1c;

% system 2 (motor to mass 2)
[S_estimate_y2c,f] = tfestimate(FRF_y2_control(:,6), FRF_y2_control(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[S_coh_y2c, ~] = mscohere(FRF_y2_control(:,6), FRF_y2_control(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_estimate_y2c,~] = tfestimate(FRF_y2_control(:,6), FRF_y2_control(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_coh_y2c, ~] = mscohere(FRF_y2_control(:,6), FRF_y2_control(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
H_estimate_y2c = (-PS_estimate_y2c)./S_estimate_y2c;
H_coh_y2c = S_coh_y2c.*PS_coh_y2c;

%%
% Create models and try to fit it to the data
m1 = 0.00021;
m2 = 0.00021;
d1 = 0.0015;
d2 = 0.0015;
k = 14;
s = tf('s');
H1_plant = (s^2*m2+d1*s+k)/(s^4*(m1*m2)+s^3*d1*(m1+m2)+s^2*k*(m1+m2));
H2_plant = (d2*s+k)/(s^4*(m1*m2)+s^3*d2*(m1+m2)+s^2*k*(m1+m2));
H21_plant = H2_plant/H1_plant;
Tdelay = 0.0018;
[num,den] = pade(Tdelay,3);
sysdelay = tf(num,den);

H21 = H21_plant*sysdelay;
H1= H1_plant*sysdelay;
H2 = H2_plant*sysdelay;

H1_model = squeeze(freqresp(H1,(f*2*pi)));
H2_model = squeeze(freqresp(H2,(f*2*pi)));
H21_model = squeeze(freqresp(H21,(f*2*pi)));

if (true)
    figure(1)
    subplot(3,1,1)
    semilogx(f,db(H_estimate_y1c), f, db(H1_model));
    grid on;
    title('Plant H1')
    ylabel('Amplitude [db]')
    legend('Measured plant','Model')
    subplot(3,1,2)
    semilogx(f,(angle(H_estimate_y1c)*180/pi),f,(angle(H1_model)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,H_coh_y1c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure(2)
    subplot(3,1,1)
    semilogx(f,db(H_estimate_y2c), f, db(H2_model));
    grid on;
    title('Plant H2')
    ylabel('Amplitude [db]')
    legend('Measured plant','Model')
    subplot(3,1,2)
    semilogx(f,(angle(H_estimate_y2c)*180/pi),f,(angle(H2_model)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,H_coh_y2c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure(3)
    bode(H21)
    title('FRF from mass 1 to mass 2')
    grid on;
end