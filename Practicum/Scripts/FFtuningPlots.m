%% FF tuning 4CM00
clear all; close all; clc;

load('FF_tune_withKfa.mat')
Kfa = Data;
load('FF_tune_withKfv.mat')
Kfv = Data;
load('FF_tune_withKfc.mat')
Kfc = Data;
load('FF_tune_withoutParameters.mat')
NoFF = Data;
clear Data;

%% Plotting
figure(1)
hold on
plot(NoFF.Time,NoFF.e,'LineWidth',0.1)
plot(Kfc.Time,Kfc.e,'LineWidth',0.1)
plot(Kfv.Time,Kfv.e,'LineWidth',0.1)
plot(Kfa.Time,Kfa.e,'LineWidth',0.1)
hold off
grid on
xlabel('Time [s]')
ylabel('Error [rad]')
xlim([0.5 8.5])
legend('No FF','Kfc tuned','Kfv tuned','Kfa tuned','Location','Best')
title('Tuning of Feed Forward Controller')

figure(2)
plot(NoFF.Time,NoFF.r)
grid on