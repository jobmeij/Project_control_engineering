%% 5CSA0 practicum session 1
clear all; close all; clc;

SampleRate = 2048;
SampleTime = 1/SampleRate;

load('data/20181012_1408.mat')
FRF_y1 = ans';
load('data/20181012_1412.mat')
FRF_y2 = ans';
load('data/20181012_1538.mat')
FRF_y1_control = ans';
load('data/20181012_1544.mat')
FRF_y2_control = ans';
load('data/data/20181012_1640.mat')
FRF_y2_FF = ans';
clear ans;
% 1. Time
% 2. enc0 [y1]
% 3. enc1 [y2]
% 4. control output (dac input [u]
% 5. error (contro input) [e]
% 6. disturbance input (inserted at control input) [d]
% 7. input/reference [r]


%% FRF data motor to encoder 1 - transfer function estimate
if (false)
    figure(1)
    hold on
    plot(FRF_y1(:,1),FRF_y1(:,2))
    plot(FRF_y1(:,1),FRF_y1(:,3))
    plot(FRF_y1(:,1),FRF_y1(:,4))
    plot(FRF_y1(:,1),FRF_y1(:,5))
    plot(FRF_y1(:,1),FRF_y1(:,6))
    plot(FRF_y1(:,1),FRF_y1(:,7))
    hold off
    grid on
    legend('Enc0','Enc1','Control out','Error','Disturbance','Input','Location','Best')
    title('Encoder 1 [y1]')
end

%% FRF data motor to encoder 2 (y2)
if (false)
    figure(2)
    hold on
    plot(FRF_y2(:,1),FRF_y2(:,2))
    plot(FRF_y2(:,1),FRF_y2(:,3))
    plot(FRF_y2(:,1),FRF_y2(:,4))
    plot(FRF_y2(:,1),FRF_y2(:,5))
    plot(FRF_y2(:,1),FRF_y2(:,6))
    plot(FRF_y2(:,1),FRF_y2(:,7))
    hold off
    grid on
    legend('Enc0','Enc1','Control out','Error','Disturbance','Input','Location','Best')
    title('Encoder 2 [y2]')
end

%% TFestimate y1
nfft = 10*(1/SampleTime); % Window length

[S_estimate_y1,f] = tfestimate(FRF_y1(:,6), FRF_y1(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[S_coh_y1, ~] = mscohere(FRF_y1(:,6), FRF_y1(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

[PS_estimate_y1,~] = tfestimate(FRF_y1(:,6), FRF_y1(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_coh_y1, ~] = mscohere(FRF_y1(:,6), FRF_y1(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

H_estimate_y1 = (-PS_estimate_y1)./S_estimate_y1;
H_coh_y1 = S_coh_y1.*PS_coh_y1;

if (false)
    figure()
    subplot(3,1,1)
    semilogx(f,db(S_estimate_y1));
    title('Sensitivity y1')
    ylabel('Amplitude [db]')
    grid on;
    subplot(3,1,2)
    semilogx(f,(angle(S_estimate_y1)*180/pi));
    ylabel('Angle [deg]')
    grid on;
    subplot(3,1,3)
    semilogx(f,S_coh_y1);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(-PS_estimate_y1));
    grid on;
    title('Process sensitivity y1')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(-PS_estimate_y1)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,PS_coh_y1);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(H_estimate_y1));
    grid on;
    title('Plant y1')
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

%% TFestimate y2
nfft = 16*(1/SampleTime); % Window length

[S_estimate_y2,f] = tfestimate(FRF_y2(:,6), FRF_y2(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[S_coh_y2, ~] = mscohere(FRF_y2(:,6), FRF_y2(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

[PS_estimate_y2,~] = tfestimate(FRF_y2(:,6), FRF_y2(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_coh_y2, ~] = mscohere(FRF_y2(:,6), FRF_y2(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

H_estimate_y2 = (-PS_estimate_y2)./S_estimate_y2;
H_coh_y2 = S_coh_y2.*PS_coh_y2;

if (false)
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
    semilogx(f,db(H_estimate_y2));
    grid on;
    title('Plant y2')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(H_estimate_y2)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,H_coh_y2);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
end

%% Tuned controller y1
if (false)
    figure()
    hold on
    plot(FRF_y1_control(:,1),FRF_y1_control(:,2))
    plot(FRF_y1_control(:,1),FRF_y1_control(:,3))
    plot(FRF_y1_control(:,1),FRF_y1_control(:,4))
    plot(FRF_y1_control(:,1),FRF_y1_control(:,5))
    plot(FRF_y1_control(:,1),FRF_y1_control(:,6))
    plot(FRF_y1_control(:,1),FRF_y1_control(:,7))
    hold off
    grid on
    legend('Enc0','Enc1','Control out','Error','Disturbance','Input','Location','Best')
    title('Encoder 1 controller [y1]')
end

nfft = 16*(1/SampleTime); % Window length

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

if (false)
    figure()
    subplot(3,1,1)
    semilogx(f,db(S_estimate_y1c));
    title('Sensitivity y1 control')
    ylabel('Amplitude [db]')
    grid on;
    subplot(3,1,2)
    semilogx(f,(angle(S_estimate_y1c)*180/pi));
    ylabel('Angle [deg]')
    grid on;
    subplot(3,1,3)
    semilogx(f,S_coh_y1c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(-PS_estimate_y1c));
    grid on;
    title('Process sensitivity y1 control')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(-PS_estimate_y1c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,PS_coh_y1c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(H_estimate_y1c));
    grid on;
    title('Plant y1 control')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(H_estimate_y1c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,H_coh_y1c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
end

%% Tuned controller y2
if (false)
    figure()
    hold on
    plot(FRF_y2_control(:,1),FRF_y2_control(:,2))
    plot(FRF_y2_control(:,1),FRF_y2_control(:,3))
    plot(FRF_y2_control(:,1),FRF_y2_control(:,4))
    plot(FRF_y2_control(:,1),FRF_y2_control(:,5))
    plot(FRF_y2_control(:,1),FRF_y2_control(:,6))
    plot(FRF_y2_control(:,1),FRF_y2_control(:,7))
    hold off
    grid on
    legend('Enc0','Enc1','Control out','Error','Disturbance','Input','Location','Best')
    title('Encoder 2 controller [y2]')
end

nfft = 16*(1/SampleTime); % Window length

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

f_y2c = f;

if (true)
    figure()
    subplot(3,1,1)
    semilogx(f,db(S_estimate_y2c));
    title('Sensitivity y2 control')
    ylabel('Amplitude [db]')
    grid on;
    subplot(3,1,2)
    semilogx(f,(angle(S_estimate_y2c)*180/pi));
    ylabel('Angle [deg]')
    grid on;
    subplot(3,1,3)
    semilogx(f,S_coh_y2c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(-PS_estimate_y2c));
    grid on;
    title('Process sensitivity y2 control')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(-PS_estimate_y2c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,PS_coh_y2c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
    
    figure()
    subplot(3,1,1)
    semilogx(f,db(H_estimate_y2c));
    grid on;
    title('Plant y2 control')
    ylabel('Amplitude [db]')
    subplot(3,1,2)
    semilogx(f,(angle(H_estimate_y2c)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    subplot(3,1,3)
    semilogx(f,H_coh_y2c);
    xlabel('Frequency [Hz]')
    ylabel('Coherence [-]')
    grid on;
end

%% Tuned feedforward controller + controller - encoder 2 (y2)
if (true)
    figure()
    hold on
    plot(FRF_y2_FF(:,1),FRF_y2_FF(:,2))
    plot(FRF_y2_FF(:,1),FRF_y2_FF(:,3))
    plot(FRF_y2_FF(:,1),FRF_y2_FF(:,4))
    plot(FRF_y2_FF(:,1),FRF_y2_FF(:,5))
    plot(FRF_y2_FF(:,1),FRF_y2_FF(:,6))
    plot(FRF_y2_FF(:,1),FRF_y2_FF(:,7))
    hold off
    grid on
    legend('Enc0','Enc1','Control out','Error','Disturbance','Input','Location','Best')
    title('Encoder 2 feedforward controller [y2]')
end


nfft = 4*(1/SampleTime); % Window length

[S_estimate_y2FF,f] = tfestimate(FRF_y2_FF(:,7), FRF_y2_FF(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[S_coh_y2FF, ~] = mscohere(FRF_y2_FF(:,7), FRF_y2_FF(:,4),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

[PS_estimate_y2FF,~] = tfestimate(FRF_y2_FF(:,7), FRF_y2_FF(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);
[PS_coh_y2FF, ~] = mscohere(FRF_y2_FF(:,7), FRF_y2_FF(:,5),...
    hann(nfft), nfft/2, nfft, 1/SampleTime);

H_estimate_y2FF = (-PS_estimate_y2FF)./S_estimate_y2FF;
H_coh_y2FF = S_coh_y2FF.*PS_coh_y2FF;

%% Determine FRF of used controller so a two point methode can be used.
% For some reason the required C files cannot be loaded when simulated from
% script. In order to simulate the model, open the model ('controller_frf.slx')
% and then set the current working directory to C:\program
% files\matlab\2017a\toolbox\shapeit and run the simulink model from
% simulink. after that the rest of the script can be used to process the
% data
if (false)
    SimTime = 1000;
    SampleRate = 2048;
    SampleTime = 1/SampleRate;
    
    %sim('controller_frf.slx');
    
    nfft = 16*(1/SampleTime); % Window length
    
    [Controller_est_SL,Controller_f] = tfestimate(Controller_sim_U_in.Data, Controller_sim_y_out.Data,...
        hann(nfft), nfft/2, nfft, 1/SampleTime);
    [Controller_coh_SL, ~] = mscohere(Controller_sim_U_in.Data, Controller_sim_y_out.Data,...
        hann(nfft), nfft/2, nfft, 1/SampleTime);
    
    Plant = (1./Controller_est_SL').*((1./S_estimate_y2c')-1);
    %S_coh_y2c
    
    if (false)
        figure()
        subplot(3,1,1)
        semilogx(Controller_f,db(Controller_est_SL));
        title('Simulated Controller')
        ylabel('Amplitude [db]')
        grid on;
        subplot(3,1,2)
        semilogx(Controller_f,(angle(Controller_est_SL)*180/pi));
        ylabel('Angle [deg]')
        grid on;
        subplot(3,1,3)
        semilogx(Controller_f,Controller_coh_SL);
        xlabel('Frequency [Hz]')
        ylabel('Coherence [-]')
        grid on;
    end
    
    if (true)
        figure()
        subplot(3,1,1)
        semilogx(f,db(Plant));
        title('Measured plant, used simulated controller')
        ylabel('Amplitude [db]')
        grid on;
        subplot(3,1,2)
        semilogx(f,(angle(Plant)*180/pi));
        ylabel('Angle [deg]')
        grid on;
        subplot(3,1,3)
        semilogx(f,S_coh_y2c);
        xlabel('Frequency [Hz]')
        ylabel('Coherence [-]')
        grid on;
    end
end

%% Transfer function estimation - process sensitivity
w = f_y2c*2*pi;                 % frequency [rad/s]
Ts = 1/2048;                    % Sample time
%EstimateOrder = 6;
EstimateOrder = 4;
x = idfrd(-PS_estimate_y2c,w,Ts);
EstimatedTf_PS = tfest(x,EstimateOrder);

boptions = bodeoptions;
boptions.FreqUnits = 'Hz';
[magPS,phasePS,woutPS] = bode(EstimatedTf_PS,boptions);

magPS = squeeze(magPS);
phasePS = squeeze(phasePS);
woutPS = squeeze(woutPS);

figure()
subplot 211
semilogx((woutPS/(2*pi)),db(magPS))
hold on
semilogx(f_y2c,db(-PS_estimate_y2c));
hold off
grid on
ylabel('Amplitude [db]')
title('Estimated vs. measured process sensitivity')
subplot 212
semilogx((woutPS/(2*pi)),phasePS)
hold on
semilogx(f_y2c,(angle(-PS_estimate_y2c)*180/pi));
hold off
grid on
ylabel('Angle [deg]')
xlabel('Frequency [Hz]')

%% TF estimation - sensitivity
w = f_y2c*2*pi;                 % frequency [rad/s]
Ts = 1/2048;                    % Sample time
%EstimateOrder = 6;
EstimateOrder = 4;
x = idfrd(S_estimate_y2c,w,Ts);
EstimatedTf_S = tfest(x,EstimateOrder);

boptions = bodeoptions;
boptions.FreqUnits = 'Hz';
[magS,phaseS,woutS] = bode(EstimatedTf_S,boptions);

magS = squeeze(magS);
phaseS = squeeze(phaseS);
woutS = squeeze(woutS);

figure()
subplot 211
semilogx((woutS/(2*pi)),db(magS))
hold on
semilogx(f_y2c,db(S_estimate_y2c));
hold off
grid on
ylabel('Amplitude [db]')
title('Estimated vs. measured sensitivity')
subplot 212
semilogx((woutS/(2*pi)),phaseS)
hold on
semilogx(f_y2c,(angle(S_estimate_y2c)*180/pi));
hold off
grid on
ylabel('Angle [deg]')
xlabel('Frequency [Hz]')

%% TF estimation - plant
TF_estimate_y2c = EstimatedTf_PS/EstimatedTf_S

% figure()
% bode(TF_estimate_y2c)
% grid on

boptions = bodeoptions;
boptions.FreqUnits = 'Hz';
[magTF,phaseTF,woutTF] = bode(TF_estimate_y2c,boptions);
magTF = squeeze(magTF);
phaseTF = squeeze(phaseTF);
woutTF = squeeze(woutTF);
foutTF = woutTF/(2*pi);

TF_complex= magTF.*cos(phaseTF*2*pi/180)+i*magTF.*sin(phaseTF*2*pi/180);



figure()
subplot 211
semilogx((woutTF/(2*pi)),db(magTF))
hold on
semilogx(f_y2c,db(H_estimate_y2c));
hold off
grid on
ylabel('Amplitude [db]')
title('Estimated vs. measured plant')
subplot 212
semilogx((woutTF/(2*pi)),phaseTF)
hold on
semilogx(f_y2c,(angle(H_estimate_y2c)*180/pi));
hold off
grid on
ylabel('Angle [deg]')
xlabel('Frequency [Hz]')

Ts = 1/2048;
DCT_TF_est_y2c = c2d(TF_estimate_y2c,Ts);

%% Controller continious to discrete
Kp = 0.5;
Kv = 0.005;
f1 = 59;
f2 = 59;
B1 = 0.03;
B2 = 0.5;

Ts = 1/2048;

CT_PD = tf([Kp],[Kv 1]);
CT_notch = tf([(1/(2*pi*f1)^2) ((2*B1)/(2*pi*f1)) 1],[(1/(2*pi*f2)^2) ((2*B2)/(2*pi*f2)) 1]);

CT_Controller = CT_PD * CT_notch;
DCT_Controller = c2d(CT_Controller,Ts);