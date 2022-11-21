%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%   4CM00 - Control Engineering     %
%   Exercise set 3                  %
%                                   %
%   Author: Marcel van Wensveen     %
%   Date: 19-09-2019                %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question 3.1

clear all, close all, clc

controllerPlot = true;
Ts = 0.0001;
simTime = 0.4;

% Plant parameters
m1 = 0.015;

m2 = 0.045;
d = 0.4;
k = 2200;
s = tf('s');
H1 = (s^2*m2+d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));
H2 = (d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));

H21 = H2/H1;


% Controller design (20Hz, 45 deg (G = 1000))
Gain = 1000;
Bw = 20; %[Hz]
PD = (1+(1/(2*pi*Bw))*s);
C = Gain*PD;

Nf1 = 70.3;     % Location of zeros
Nf2 = 35.2;     % Location of poles
Nbeta1 = 0.1;     
Nbeta2 = 1;
Notch = ((1/(2*pi*Nf1)^2)*s^2+(2*Nbeta1/(2*pi*Nf1))*s+1)/((1/(2*pi*Nf2)^2)*s^2+(2*Nbeta2/(2*pi*Nf2))*s+1);

boptions = bodeoptions;
boptions.FreqUnits = 'Hz';

if controllerPlot 
%     figure(1)
%     bode(H1)
%     title('H1')
%     figure(2)
%     bode(C)
%     title('C')
    figure(3)
    bode(H1*C*Notch)
    %title('H*C')
    title('Open loop transfer function - F to x1')
    grid on
    figure(4)
    bode(1/(1+H1*C))
    title('S')
    figure(5)
    nyquist(H1*C)
    title('H*C')
    figure(20)
    step((H1*C)/(1+H1*C))
    title('Step response')
end
    
sim('exercise_3_1.slx')
CL_Output_bw20 = y_out;
error_Output_bw20 = error_out;

% Controller 2 (Bw variation) 10Hz, 45 deg
Gain = 200;
Bw = 10; %[Hz]
PD = (1+(1/(2*pi*Bw))*s);
C = Gain*PD;
sim('exercise_3_1.slx')
CL_Output_bw10 = y_out;
error_Output_bw10 = error_out;

% Controller 3 (Bw variation) 5Hz, 45 deg
Gain = 48;
Bw = 5; %[Hz]
PD = (1+(1/(2*pi*Bw))*s);
C = Gain*PD;
sim('exercise_3_1.slx')
CL_Output_bw5 = y_out;
error_Output_bw5 = error_out;

% Plot Bw variation 
figure(6)
plot(Step_In)
hold on
plot(CL_Output_bw5)
plot(CL_Output_bw10)
plot(CL_Output_bw20)
hold off
legend('Step input','Bw = 5 Hz','Bw = 10 Hz','Bw = 20 Hz')
title('Step response with different bandwidth, phase margin = 35 [deg]')
xlabel('Time [sec]')
ylabel('Amplitude [-]')
figure(7)
plot(Step_In)
hold on
plot(error_Output_bw5)
plot(error_Output_bw10)
plot(error_Output_bw20)
hold off
legend('Step input','Bw = 5 Hz','Bw = 10 Hz','Bw = 20 Hz')
title('Error to step response with different bandwidth, phase margin = 35 [deg]')
xlabel('Time [sec]')
ylabel('Amplitude [-]')

% Controller 4 (PM variation) 20hz 45 deg
Gain = 950;
Bw = 21; %[Hz]
PD = (1+(1/(2*pi*Bw))*s);
C = Gain*PD;
sim('exercise_3_1.slx')
CL_Output_PM45 = y_out;
error_Output_PM45 = error_out;

% Controller 5 (PM variation) 20Hz, 35 deg
Gain = 1100;
Bw = 30; %[Hz]
PD = (1+(1/(2*pi*Bw))*s);
C = Gain*PD;
sim('exercise_3_1.slx')
CL_Output_PM35 = y_out;
error_Output_PM35 = error_out;

% Controller 6 (PM variation) 20Hz, 55 deg
Gain = 790;
Bw = 15; %[Hz]
PD = (1+(1/(2*pi*Bw))*s);
C = Gain*PD;
sim('exercise_3_1.slx')
CL_Output_PM55 = y_out;
error_Output_PM55 = error_out;

% Plot Bw variation 
figure(8)
plot(Step_In)
hold on
plot(CL_Output_PM35)
plot(CL_Output_PM45)
plot(CL_Output_PM55)
hold off
legend('Step input','PM = 35 [deg]','PM = 45 [deg]','PM = 55 [deg]')
title('Step response with different phase margin, bandwidth = 20 Hz')
xlabel('Time [sec]')
ylabel('Amplitude [-]')
figure(9)
plot(Step_In)
hold on
plot(error_Output_PM35)
plot(error_Output_PM45)
plot(error_Output_PM55)
hold off
legend('Step input','PM = 35 [deg]','PM = 45 [deg]','PM = 55 [deg]')
title('Error to step response with different phase margin, bandwidth = 20 Hz')
xlabel('Time [sec]')
ylabel('Amplitude [-]')

%% test controller to H2 (question 3.1)

H2 = (d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));

% Controller design (20Hz, 45 deg (G = 1000))
Gain = 1000;
Bw = 20; %[Hz]
PD = (1+(1/(2*pi*Bw))*s);
C = Gain*PD;

figure(1)
bode(H2)
title('H2')
figure(2)
bode(C)
title('C')
figure(3)
bode(H2*C)
title('H2*C')
figure(4)
bode(1/(1+H2*C))
title('S2')
figure(5)
nyquist(H2*C)
title('H2*C')
figure(20)
step((H2*C)/(1+H2*C))
title('Step response')


%% Question 3.2 - Inverted pendulum

clear all, close all, clc

M = 0.05; % [kg]
m = 0.04; % [kg]
b = 0.02; % [s/m]
g = 9.81; % [m/s^2]
L = 0.15; % [m]

s = tf('s');

H = s/(((-(4*M+m)*L*s^3)/6)-((2*L*b*s^2)/3)+(g*(M+m)*s)+g*b);

%% part a
figure(1)
bode(H);
grid on;

H_zeros = zero(H);
H_poles  = pole(H);
disp('Open loop poles of H located at: ');
disp(H_poles);
disp('Open loop zero of H located at: ');
disp(H_zeros);

%% part b
C1 = 1;
C2 = -0.55;

figure(2)
nyquist(C1*H)
title('C1*H (Gain = 1)')
figure(3)
nyquist(C2*H)
title('C2*H (Gain = -0.55)')

figure(4)
step((C1*H)/(1+C1*H))
title('Step response of C1*H (Gain = 1)')
figure(5)
step((C2*H)/(1+C2*H))
title('Step response of C2*H (Gain = -0.55)')

%% part c

Gain = -4;
Bw = 7; %[Hz] = 62.8 rad
D = (1+(1/(2*pi*Bw))*s);
I = (s+(2*pi*Bw/5))/s;
LP = 1/((s/(2*pi*Bw*5))+1);
Lead = (((s/(2*pi*Bw/3))+1)/((s/(2*pi*Bw*3))+1));
Lag = (((s/(2*pi*Bw*3))+1)/((s/(2*pi*Bw/3))+1));
b1 = 1;
b2 = 10;
notchLead = (((s^2/(2*pi*Bw/3)^2)+((2*b1*s)/(2*pi*Bw/3))+1)/((s^2/(2*pi*Bw*3)^2)+((2*b2*s)/(2*pi*Bw*3))+1));
notchLag = (((s^2/(2*pi*Bw*3)^2)+((2*b1*s)/(2*pi*Bw*3))+1)/((s^2/(2*pi*Bw/3)^2)+((2*b2*s)/(2*pi*Bw/3))+1));
notchFreq = 20.9; % [rad]
b1_notch = 1;
b2_notch = 1;
notch = (((s^2/notchFreq^2)+((2*b1_notch*s)/notchFreq)+1)/((s^2/notchFreq^2)+((2*b2_notch*s)/notchFreq)+1));
C = Gain*Lead*D*I*I*LP;

n = 1000; %// Define number of points on circle
theta = linspace(0, 2*pi, n);
x1 = 0.5*cos(theta);
y1 = 0.5*sin(theta);
x2 = cos(theta);
y2 = sin(theta);
figure(6)
nyquist(C*H)
hold on
plot(x1-1, y1, 'r', x2, y2, 'b'); %// Unit circle
hold off
axis([-5 5 -5 5]);
title('Nyquist of C*H')
figure(7)
step((C*H)/(1+C*H))
grid on
title('CL step response')
figure(8)
bode(C*H)
title('Bode of C*H')
grid on;
figure(9)
bode(1/(1+C*H))
title('Bode of S')
grid on;

%% Simulink model
Disturbance = 1;
stepTime = 1;
simTime = 2;
Ts = 0.001;

sim('exercise_3_2.slx')

figure(10)
plot(y_out)
hold on;
plot(x_in)
hold off;
xlabel('Time [s]')
ylabel('Amplitude [-]')
legend('Output position','Input reference')
title('Step response with constant disturbance')

%% Exercise 3.3

s = tf('s');
H1 = 1000*((s+6)/(s^3+20*s^2+5000*s));
H2 = 1000*((s-6)/(s^3+20*s^2+5000*s));

Gain = 20;
Bw = 30; %[Hz] 30 = 188.4 rad
D = (1+(1/(2*pi*Bw))*s);
I = (s+(2*pi*Bw/5))/s;
LP = 1/((s/(2*pi*Bw*5))+1);
Lead = (((s/(2*pi*Bw/3))+1)/((s/(2*pi*Bw*3))+1));
Lag = (((s/(2*pi*Bw*3))+1)/((s/(2*pi*Bw/3))+1));

C1 = Gain*D;

n = 1000; %// Define number of points on circle
theta = linspace(0, 2*pi, n);
x1 = 0.67*cos(theta);
y1 = 0.67*sin(theta);
x2 = cos(theta);
y2 = sin(theta);
figure(6)
nyquist(C1*H1)
hold on
plot(x1-1, y1, 'r', x2, y2, 'b'); %// Unit circle
hold off
axis([-4 4 -4 4]);
title('Nyquist of C1*H1')
figure(7)
step((C1*H1)/(1+C1*H1))
grid on
title('C1*H1 CL step response')
figure(8)
bode(C1*H1)
title('Bode of C1*H1')
grid on;
figure(9)
bode(1/(1+C1*H1))
title('Bode of S1')
grid on;


figure(10)
bode(C1*H2, C1*H2)
title('Bode of C1*H2')
figure(11)
nyquist(C1*H2)
hold on
plot(x1-1, y1, 'r', x2, y2, 'b'); %// Unit circle
hold off
axis([-4 4 -4 4]);
title('Nyquist of C1*H2')
figure(12)
step((C1*H2)/(1+C1*H2))
grid on
title('C1*H2 CL step response')

%% 3.3 for H2
s = tf('s');
H2 = 1000*((s-6)/(s^3+20*s^2+5000*s));

Gain = 600000;
Bw = 150; %[Hz] 30 = 188.4 rad
D = (1+(1/(2*pi*Bw))*s);
I = (s+(2*pi*Bw/5))/s;
LP = 1/((s/(2*pi*Bw*50))+1);
Lead = (((s/(2*pi*Bw/3))+1)/((s/(2*pi*Bw*3))+1));
Lag = (((s/(2*pi*Bw*3))+1)/((s/(2*pi*Bw/3))+1));
b1 = 1;
b2 = 10;
notchLead = (((s^2/(2*pi*Bw/3)^2)+((2*b1*s)/(2*pi*Bw/3))+1)/((s^2/(2*pi*Bw*3)^2)+((2*b2*s)/(2*pi*Bw*3))+1));
notchLag = (((s^2/(2*pi*Bw*3)^2)+((2*b1*s)/(2*pi*Bw*3))+1)/((s^2/(2*pi*Bw/3)^2)+((2*b2*s)/(2*pi*Bw/3))+1));
notchFreq = 20.9; % [rad]
b1_notch = 1;
b2_notch = 1;
notch = (((s^2/notchFreq^2)+((2*b1_notch*s)/notchFreq)+1)/((s^2/notchFreq^2)+((2*b2_notch*s)/notchFreq)+1));

C2 = Gain*D*notchLag*LP;

n = 1000; %// Define number of points on circle
theta = linspace(0, 2*pi, n);
x1 = 0.5*cos(theta);
y1 = 0.5*sin(theta);
x2 = cos(theta);
y2 = sin(theta);
figure(6)
nyquist(C2*H2)
hold on
plot(x1-1, y1, 'r', x2, y2, 'b'); %// Unit circle
hold off
axis([-4 4 -4 4]);
title('Nyquist of C2*H2')
figure(7)
step((C2*H2)/(1+C2*H2))
grid on
title('C2*H2 CL step response')
figure(8)
bode(C2*H2)
title('Bode of C2*H2')
grid on;
figure(9)
bode(1/(1+C2*H2))
title('Bode of S2')
grid on;
