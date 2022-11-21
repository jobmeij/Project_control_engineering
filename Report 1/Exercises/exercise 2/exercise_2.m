%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%   4CM00 - Control Engineering     %
%   Exercise set 2                  %
%                                   %
%   Author: Marcel van Wensveen     %
%   Date: 19-09-2019                %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question 2.1

clear all, close all, clc

f1 = 40; % [Hz]
f2 = 40; % [Hz]
b1 = 0.01;
b2 = 1;
Ts = 0.001;
Fs = 1/Ts;
TimeVector = [0:Ts:30]';
num = [(1/((2*pi*f1)^2)) ((2*b1)/(2*pi*f1)) 1];
den = [(1/((2*pi*f2)^2)) ((2*b2)/(2*pi*f2)) 1];
%Notch = tf(num, den);
%bode(Notch);

sim('exercise_2_1.slx');

figure(1)
hold on
plot(TimeVector,BeforeFiltering)
% figure()
plot(TimeVector,AfterFiltering)
title('Random noise before and after notch filter')
legend('Before filtering','After filtering')
xlabel('Time [s]')
ylabel('Amplitude [-]')
hold off


fftBeforeFiltering = fft(BeforeFiltering);
spectrumBeforeFilter = abs(fftBeforeFiltering/length(fftBeforeFiltering));
fftAfterFiltering = fft(AfterFiltering);
spectrumAfterFilter = abs(fftAfterFiltering/length(fftAfterFiltering));
freqVector = linspace(0, Fs/2, (length(fftBeforeFiltering)/2)+1);

CrossPower = spectrumAfterFilter./spectrumBeforeFilter;

figure(2)
hold on
plot(freqVector, db(spectrumBeforeFilter(1:round((length(fftBeforeFiltering)/2)))));
plot(freqVector, db(spectrumAfterFilter(1:round((length(fftAfterFiltering)/2)))));
legend('Before filtering','After filtering')
xlabel('Frequency [Hz]')
ylabel('Amplitude [db]')
title('FFT of signals')
hold off

figure(3)
plot(freqVector, db(CrossPower(1:round((length(fftAfterFiltering)/2)))));
xlabel('Frequency [Hz]')
ylabel('Amplitude [db]')
title('CPSD of notch filter')

%% Question 2.2

clear all, close all, clc

Bw = 10*2*pi; % [rad]

p_gain = 2.7887e+03;
d_gain = 4.4384e+01;
m = 1;

Ts = 0.001;
nfft = 30*(1/Ts); % Window length

H = tf(1,[m 0 0]);
C = tf([d_gain p_gain],1);
%bode(H*C);

SimulationTime = 120; %[sec]
sim('exercise_2_2.slx')
TimeVector = [0:0.001:120]';

figure(1)
hold on
plot(measured_d{1}.Values)
plot(measured_u{1}.Values)
hold off
title('Time domain plots of d and u')
legend('d (input noise)', 'u (input to plant)')

fft_d = fft(measured_d{1}.Values.Data);
fft_u = fft(measured_u{1}.Values.Data);
freqVector = linspace(0,500,(length(fft_u)/2)+1);

figure(2)
plot(freqVector,db(fft_d(1:round(length(fft_u)/2))));
title('u (input to plant)')
figure(3)
plot(freqVector,db(fft_u(1:round(length(fft_d)/2))));
title('d (input noise)')

[Txy,f] = tfestimate(measured_d{1}.Values.Data,...
    measured_u{1}.Values.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);
[Cxy, ~] = mscohere(measured_d{1}.Values.Data,...
    measured_u{1}.Values.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);

Data = freqresp(C, f, 'Hz');
Cfft = Data(:)';

L_estimate = ((1./Txy)-1);
H_estimate = (1./(Cfft))'.*((1./Txy)-1);


figure(4)
subplot(3,1,1)
semilogx(f, db(Txy))
title('Sensitivity measurement')
subplot(3,1,2)
semilogx(f, angle(Txy)*180/pi)
subplot(3,1,3)
semilogx(f, Cxy)

figure(5)
subplot(2,1,1)
semilogx(f, db(Cfft))
title('Controller frequency behaviour')
subplot(2,1,2)
semilogx(f, angle(Cfft)*180/pi)

figure(6)
subplot(3,1,1)
semilogx(f, db(L_estimate))
grid on;
title('L (C*H) estimate')
subplot(3,1,2)
semilogx(f, angle(L_estimate)*180/pi)
grid on;
subplot(3,1,3)
semilogx(f, Cxy)
grid on;

figure(7)
subplot(3,1,1)
semilogx(f, db(H_estimate))
grid on;
title('H estimate')
subplot(3,1,2)
semilogx(f, angle(H_estimate)*180/pi)
grid on;
subplot(3,1,3)
semilogx(f, Cxy)
grid on;

%% Question 2.3

clear all, close all, clc

Ts = 1/1000;
SimTime = 15;
nfft = 4*(1/Ts); % Window length

sim('frf_ex3.mdl')

[S_estimate,f] = tfestimate(d.Data, u.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);
[S_coh, ~] = mscohere(d.Data, u.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);

[PS_estimate,~] = tfestimate(d.Data, e.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);
[PS_coh, ~] = mscohere(d.Data, e.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);

H_estimate = (-PS_estimate)./S_estimate;
H_coh = S_coh.*PS_coh;

figure(1)
subplot(3,1,1)
semilogx(f,db(S_estimate));
title('Sensitivity')
ylabel('Amplitude [db]')
grid on;
subplot(3,1,2)
semilogx(f,(angle(S_estimate)*180/pi));
ylabel('Angle [deg]')
grid on;
subplot(3,1,3)
semilogx(f,S_coh);
xlabel('Frequency [Hz]')
ylabel('Coherence [-]')
grid on;

figure(2)
subplot(3,1,1)
semilogx(f,db(-PS_estimate));
grid on;
title('Process sensitivity')
ylabel('Amplitude [db]')
subplot(3,1,2)
semilogx(f,(angle(-PS_estimate)*180/pi));
grid on;
ylabel('Angle [deg]')
subplot(3,1,3)
semilogx(f,PS_coh);
xlabel('Frequency [Hz]')
ylabel('Coherence [-]')
grid on;

figure(3)
subplot(3,1,1)
semilogx(f,db(H_estimate));
grid on;
title('Plant')
ylabel('Amplitude [db]')
subplot(3,1,2)
semilogx(f,(angle(H_estimate)*180/pi));
grid on;
ylabel('Angle [deg]')
subplot(3,1,3)
semilogx(f,H_coh);
xlabel('Frequency [Hz]')
ylabel('Coherence [-]')
grid on;

H_estimate4 = H_estimate;
H_coh4 = H_coh;
f4 = f;

%% change nfft value

nfft = 15*(1/Ts); % Window length

[S_estimate,f] = tfestimate(d.Data, u.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);
[S_coh, ~] = mscohere(d.Data, u.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);

[PS_estimate,~] = tfestimate(d.Data, e.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);
[PS_coh, ~] = mscohere(d.Data, e.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);

H_estimate15 = (-PS_estimate)./S_estimate;
H_coh15 = S_coh.*PS_coh;
f15 = f;

nfft = 1*(1/Ts); % Window length

[S_estimate,f] = tfestimate(d.Data, u.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);
[S_coh, ~] = mscohere(d.Data, u.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);

[PS_estimate,~] = tfestimate(d.Data, e.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);
[PS_coh, ~] = mscohere(d.Data, e.Data,...
    hann(nfft), nfft/2, nfft, 1/Ts);

f1 = f;
H_estimate1 = (-PS_estimate)./S_estimate;
H_coh1 = S_coh.*PS_coh;

% make sure the vectors are of the same length
H_1 = interp1(f1,H_estimate1, f15);
Coh_1 = interp1(f1,H_coh1, f15);
H_4 = interp1(f4,H_estimate4, f15);
Coh_4 = interp1(f4,H_coh4, f15);

corr_Index = find(f15 == f1(2));
H_1(1:corr_Index) = nan(corr_Index,1);
Coh_1(1:corr_Index) = nan(corr_Index,1);
corr_Index = find(f15 >= f4(2), 1);
H_4(1:corr_Index) = nan(corr_Index,1);
Coh_4(1:corr_Index) = nan(corr_Index,1);

% Create delay function to plot in same plot
s = tf('s');
sys = exp(-0.001*s);  % Sampling delay of 0.001 second  
sysx = pade(sys,3)
delayFRF = freqresp(sysx,f15, 'Hz');
H_delay = squeeze(delayFRF);

% Plot bodes
figure(4)
subplot(3,1,1)
semilogx(f15,db(H_estimate15),f15,db(H_4),f15,db(H_1), f15, nan(length(f15),1));
grid on;
title('Estimated plant with different values for nfft and estimated sampling delay')
ylabel('Amplitude [db]')
legend('nfft = 15 sec','nfft = 4 sec','nfft = 1 sec', 'Sampling delay' )
subplot(3,1,2)
semilogx(f15,(angle(H_estimate15)*180/pi),f15,(angle(H_4)*180/pi),f15,(angle(H_1)*180/pi),f15,180+(angle(H_delay)*180/pi));
grid on;
ylabel('Angle [deg]')
subplot(3,1,3)
semilogx(f15,H_coh15,f15,Coh_4,f15,Coh_1);
xlabel('Frequency [Hz]')
ylabel('Coherence [-]')
grid on;

% Time delay simulation

s = tf('s');
sys = exp(-0.1*s);    
sysx = pade(sys,3)
delayFRF = freqresp(sysx,f15, 'Hz');
H_delay = squeeze(delayFRF);

