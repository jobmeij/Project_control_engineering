%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%   4CM00 - Control Engineering     %
%   Exercise set 1                  %
%                                   %
%   Author: Marcel van Wensveen     %
%   Date: 14-09-2019                %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Question 1.1.g)

clear all, close all, clc

s = tf([1 0],[1]);

P = (s+2)^2/(s^2*(s+2)*(s^2 + s+25));
%bode(P)
omegaVector = logspace(-1,1,1e3); 

freqResponseP1 = freqresp(P,omegaVector, 'Hz');
freqResponseP = freqResponseP1(:);

% Plot data
figure()
subplot 211
semilogx(omegaVector, db(freqResponseP))
xlabel('Frequency [Hz]');
ylabel('Magnitude [db]');
grid on;
subplot 212
semilogx(omegaVector, (180/pi)*angle(freqResponseP))
xlabel('Frequency [Hz]');
ylabel('Phase [Hz]');
grid on

%% Question 1.3

clear all, close all, clc

load('frfdata.mat')

HzVector = logspace(-1,3,1e3); 

% Try to fit a TF to the data
[TF1.num, TF1.den] = frfit(H1, hz, [6,4,0], 1);
[TF2.num, TF2.den] = frfit(H2, hz, [8,2,0], 1);

TF1.TF = tf(TF1.num, TF1.den);
TF2.TF = tf(TF2.num, TF2.den);

TF1.freqResponseRaw = freqresp(TF1.TF,HzVector, 'Hz');
TF1.freqResponse = TF1.freqResponseRaw(:);
%delete(TF1.freqResponseRaw);

TF2.freqResponseRaw = freqresp(TF2.TF,HzVector, 'Hz');
TF2.freqResponse = TF2.freqResponseRaw(:);
%delete(TF2.freqResponseRaw);

% Plot data
figure()
subplot 221
semilogx(hz, db(H1), HzVector, db(TF1.freqResponse))
%xlabel('Frequency [Hz]');
ylabel('Magnitude [db]');
title('H1')
grid on;
subplot 223
semilogx(hz, (180/pi)*angle(H1),  HzVector, (180/pi)*angle(TF1.freqResponse))
xlabel('Frequency [Hz]');
ylabel('Phase [Hz]');
grid on
legend('Measured data','Estimated bode')

subplot 222
semilogx(hz, db(H2), HzVector, db(TF2.freqResponse))
%xlabel('Frequency [Hz]');
ylabel('Magnitude [db]');
title('H2')
grid on;
subplot 224
semilogx(hz, (180/pi)*angle(H2), HzVector, (180/pi)*angle(TF2.freqResponse))
xlabel('Frequency [Hz]');
ylabel('Phase [Hz]');
grid on
legend('Measured data','Estimated bode')



% Calculate pole locations of the found TF
PolesH1 = sign(real(roots(TF1.den))).*abs(roots(TF1.den)); % [rad]

figure()
subplot 211
semilogx(hz,db(H1))
hold on
semilogx(hz,db(H2))
hold off
grid on
ylabel('Gain [dB]')
title('System H1 and H2 from frfdata.mat')
legend('H1','H2')
xlim([hz(1,1) hz(end,1)])

subplot 212
semilogx(hz,(angle(H1))*(180/pi))
hold on
semilogx(hz,(angle(H2))*(180/pi))
hold off
grid on
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
xlim([hz(1,1) hz(end,1)])

%% Question 1.5

clear all, close all, clc

m1 = 0.015;
m2 = 0.045;
d = 0.4;
k = 2200;

%s = tf([1 0], [1]);
s = tf('s')

H1 = (s^2*m2+d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));

H2 = (d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));
bode(H1,H2)
