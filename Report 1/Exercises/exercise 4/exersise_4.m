%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%   4CM00 - Control Engineering     %
%   Exercise set 3                  %
%                                   %
%   Author: Marcel van Wensveen     %
%   Date: 19-09-2019                %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question 4.1

clear all, close all, clc


%% Question 4.2
% everyting is done with shapeit
% final controller is stored in 
% "exer_4_2_shapeit_controller.mat"

clear all, close all, clc

load('exer_4_2_shapeit_controller.mat');

H = shapeit_data.P.sys;
C = shapeit_data.C_tf;

options = bodeoptions;
options.FreqUnits = 'Hz';

figure(1)
bode(H, options)
title('Bode of plant (H)')
grid on;
figure(2)
bode(C, options)
title('Bode of controller (C)')
grid on;
figure(3)
bode(C*H, options)
title('Bode of open loop (C*H)')
grid on;
figure(4)
step(C*H/(1+C*H))
figure(5)
bode(H/(1+C*H))

figure(6)
lsim((H/(1+H*C)))
