%% 4CM00 exercise 3.2
clear all; close all; clc;

% Plant parameters
M = 0.05;
m = 0.04;
b = 0.02;
g = 9.81;
L = 0.15;

% Plant TF
H = tf([1 0],[-(1/6)*(4*M+m)*L -(2/3)*L*b g*(M+m) g*b])

% Determine poles and zeros
zeros = zero(H)
poles = pole(H)

% Controllers
C1 = 1;
C2 = -0.55;
C3 = tf([1],[1 0 0]) * tf([1/(-12.076110925) 1],[1]) * tf([1/(0.2222) 1],[1]) * tf([12.1873 1],[1])

%% Working controller
C4 = -1 * tf([1],[1 0 0]) * tf([1 0.2222],[1])

% Notch
s = tf('s')
Nf1 = 38.4/(2*pi);     % Location of zeros
Nf2 = 38.4/(2*pi);     % Location of poles
Nbeta1 = 1;     
Nbeta2 = 1;
Notch = ((1/(2*pi*Nf1)^2)*s^2+(2*Nbeta1/(2*pi*Nf1))*s+1)/((1/(2*pi*Nf2)^2)*s^2+(2*Nbeta2/(2*pi*Nf2))*s+1);
bode(Notch);
grid on

% Resonance 2nd order
PeakFreq1_rads = 1;
Beta1 = 0.01;
Resonance1 = tf([PeakFreq1_rads^2],[1 2*Beta1*PeakFreq1_rads PeakFreq1_rads^2])
%
PeakFreq2_rads = 1;
Beta2 = 10;
Resonance2 = tf([PeakFreq2_rads^2],[1 2*Beta2*PeakFreq2_rads PeakFreq2_rads^2])

Gain_dB = -16.2 + 3.7;

zero_freq = 10;
zero_Control = tf([1 zero_freq],[1]) * tf([1 zero_freq],[1])

C4 = C4 * 10^(Gain_dB/20) * zero_Control * Notch %* Resonance1;

bode(C4*H)
grid on

T4 = (C4*H)/(1+C4*H);
step(T4)
grid on


%% Closed loop TF 
T1 = (C1*H)/(1+C1*H);
T2 = (C2*H)/(1+C2*H);
T3 = (C3*H)/(1+C3*H);


% Plotting things
boptions = bodeoptions;
boptions.FreqUnits = 'Hz';

figure(1)
bode(H,boptions)
grid on

figure(2)
nyquist(C1*H)
hold on
nyquist(C2*H)
grid on
legend('C1 H','C2 H')

figure(3)
step(T1)
hold on
step(T2)
grid on
legend('Step T1','Step T2')

figure(4)
step(T3)
grid on

