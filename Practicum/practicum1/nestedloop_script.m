%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
%   4CM00 - Control Engineering      %
% Practicum - Nested loop simulation %
%                                    %
%   Author: Marcel van Wensveen      %
%       &   Job Meijer               %
%   Date: 18-10-2019                 %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create models

%clear all; close all; clc;

f = linspace(0,1024,16385)';

% Create models and try to fit it to the data
Tdelay = 0.0018;
m1 = 0.00021;
m2 = 0.00021;
d = 0.0015;
k = 14;
s = tf('s');
H1_plant = (s^2*m2+d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));
H2_plant = (d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));
H21_plant = H2_plant/H1_plant;

[num,den] = pade(Tdelay,3);
sysdelay = tf(num,den);

H21 = H21_plant*sysdelay;
H1= H1_plant*sysdelay;
H2 = H2_plant*sysdelay;

H1_model = squeeze(freqresp(H1,(f*2*pi)));
H2_model = squeeze(freqresp(H2,(f*2*pi)));
H21_model = squeeze(freqresp(H21,(f*2*pi)));

%% Simulink model with nested loop 

% Controller for mass 1

% gain
H1_gain = 130; %130
% notch
Nf1 = 58.1;     % Location of zeros [hz]
Nf2 = 41.1;     % Location of poles [hz]
Nbeta1 = 0.02;     
Nbeta2 = 0.013;
H1_notch = ((1/(2*pi*Nf1)^2)*s^2+(2*Nbeta1/(2*pi*Nf1))*s+1)/((1/(2*pi*Nf2)^2)*s^2+(2*Nbeta2/(2*pi*Nf2))*s+1);
% PD
P = 0.05;
D = 0.001;
H1_pd = P+D*s;
% Lowpass 1st order
f_LP = 500;
H1_LP = 1/((1/(2*pi*f_LP))*s+1);
% combine blocks
C1 = H1_gain*H1_notch*H1_pd*H1_LP;
C1_resp = squeeze(freqresp(C1,f*2*pi));

if (true)
    figure(4)
    subplot(2,1,1)
    semilogx(f,db(C1_resp),f,db(H1_model),f,db(C1_resp.*H1_model));
    grid on;
    title('System from motor to mass 1')
    ylabel('Amplitude [db]')
    legend('C1','H1','C1*H1')
    subplot(2,1,2)
    semilogx(f,(angle(C1_resp)*180/pi),f,(angle(H1_model)*180/pi),f,(angle(C1_resp.*H1_model)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    xlabel('Frequency [hz]')
end

% Controller for H21 (from mass 1 to mass 2)

% gain
H21_gain = 15; %15
% Integrator
I_I = 1;
H21_I = (1+(I_I*2*pi/s));
% Notch
Nf21 = 40.75;     % Location of zeros [hz]
Nf22 = 150;     % Location of poles [hz]
Nbeta21 = 0.01;     
Nbeta22 = 0.1;
H21_notch = ((1/(2*pi*Nf21)^2)*s^2+(2*Nbeta21/(2*pi*Nf21))*s+1)/((1/(2*pi*Nf22)^2)*s^2+(2*Nbeta22/(2*pi*Nf22))*s+1);
% Lowpass 1
f1_LP21 = 10;
H21_LP1 = 1/((1/(2*pi*f1_LP21))*s+1);
% Lowpass 2
f2_LP21 = 0.5;
H21_LP2 = 1/((1/(2*pi*f2_LP21))*s+1);
% combine blocks
C21 = H21_gain*H21_notch*H21_LP1*H21_LP2*H21_I;
C21_resp = squeeze(freqresp(C21,f*2*pi));

if (true)
    figure(5)
    subplot(2,1,1)
    semilogx(f,db(C21_resp),f,db(H21_model),f,db(C21_resp.*H21_model));
    grid on;
    title('System from mass 1 to mass 2')
    ylabel('Amplitude [db]')
    legend('C21','H21','C21*H21')
    subplot(2,1,2)
    semilogx(f,(angle(C21_resp)*180/pi),f,(angle(H21_model)*180/pi),f,(angle(C21_resp.*H21_model)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    xlabel('Frequency [hz]')
end

%% FRF from simulation to check if everyting is correct
% make sure the correct noise inputs are set in the model

% Init data
SimTime = 5; % [sec]
SampleRate = 2048;
SampleTime = 1/SampleRate;
nfft_sim = 5*(1/SampleTime/10); % Window length

disp('Starting simulation')
sim('nestedloop_simulation_jobmei');
disp('Simulation finished')


% system 1 (motor to mass 1)
[Sim_S_estimate_y1c,f_sim] = tfestimate(sim_d.data, sim_u.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
[Sim_S_coh_y1c, ~] = mscohere(sim_d.data, sim_u.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
[Sim_PS_estimate_y1c,~] = tfestimate(sim_d.data, sim_e.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
[Sim_PS_coh_y1c, ~] = mscohere(sim_d.data, sim_e.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
Sim_H1_estimate_y1c = (-Sim_PS_estimate_y1c)./Sim_S_estimate_y1c;
Sim_H1_coh_y1c = Sim_S_coh_y1c.*Sim_PS_coh_y1c;

if (false) % FRF plots for system 1 
    figure(6)
    subplot(3,1,1)
    semilogx(f_sim,db(Sim_S_estimate_y1c));
    title('Sensitivity y1')
    ylabel('Amplitude [db]')
    grid on;
    subplot(3,1,2)
    semilogx(f_sim,(angle(Sim_S_estimate_y1c)*180/pi));
    ylabel('Angle [deg]')
    grid on;
    subplot(3,1,3)
    semilogx(f_sim,Sim_S_coh_y1c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure(7)
    subplot(3,1,1)
    semilogx(f_sim,db(-Sim_PS_estimate_y1c));
    grid on;
    title('Process sensitivity y1')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f_sim,(angle(-Sim_PS_estimate_y1c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f_sim,Sim_PS_coh_y1c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure(8)
    subplot(3,1,1)
    semilogx(f_sim,db(Sim_H1_estimate_y1c));
    grid on;
    title('Simulated plant H1')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f_sim,(angle(Sim_H1_estimate_y1c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f_sim,Sim_H1_coh_y1c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
end


% system 21 (motor to mass 2, with nested loop)
[Sim_S_estimate_y21c,~] = tfestimate(sim21_d.data, sim21_u.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
[Sim_S_coh_y21c, ~] = mscohere(sim21_d.data, sim21_u.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
[Sim_PS_estimate_y21c,~] = tfestimate(sim21_d.data, sim21_e.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
[Sim_PS_coh_y21c, ~] = mscohere(sim21_d.data, sim21_e.data,...
    hann(nfft_sim), nfft_sim/2, nfft_sim, 1/SampleTime);
Sim_H_estimate_y21c = (-Sim_PS_estimate_y21c)./Sim_S_estimate_y21c;
Sim_H_coh_y21c = Sim_S_coh_y21c.*Sim_PS_coh_y21c;

if (false) % FRF plots for system 21
    figure(9)
    subplot(3,1,1)
    semilogx(f_sim,db(Sim_S_estimate_y21c));
    title('Sensitivity y21')
    ylabel('Amplitude [db]')
    grid on;
    subplot(3,1,2)
    semilogx(f_sim,(angle(Sim_S_estimate_y21c)*180/pi));
    ylabel('Angle [deg]')
    grid on;
    subplot(3,1,3)
    semilogx(f_sim,Sim_S_coh_y21c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure(10)
    subplot(3,1,1)
    semilogx(f_sim,db(-Sim_PS_estimate_y21c));
    grid on;
    title('Process sensitivity y21')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f_sim,(angle(-Sim_PS_estimate_y21c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f_sim,Sim_PS_coh_y21c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure(11)
    subplot(3,1,1)
    semilogx(f_sim,db(Sim_H_estimate_y21c));
    grid on;
    title('Simulated plant H21')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f_sim,(angle(Sim_H_estimate_y21c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f_sim,Sim_H_coh_y21c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
end


