%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%   4CM00 - Control Engineering                 %
%   Practicum - Single controller simulation    %
%                                               %
%   Author: Marcel van Wensveen                 %
%       &   Job Meijer                          %
%   Date:   17-10-2019                          %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Create models and try to fit it to the data
SampleRate = 2048;
SampleTime = 1/SampleRate;

f = linspace(0,1024,16385)';

m1 = 0.00021;
m2 = 0.00021;
d = 0.0015;
k = 14;
s = tf('s');
H1_plant = (s^2*m2+d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));
H2_plant = (d*s+k)/(s^4*(m1*m2)+s^3*(d*m1+d*m2)+s^2*(m1*k+m2*k));
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

if (false)
    figure(1)
    subplot(3,1,1)
    semilogx(f,db(H_estimate_y1c), f, db(H1_model));
    grid on;
    title('Measured plant Y1')
    ylabel('Amplitude [db]')
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
    title('Measured plant Y2')
    ylabel('Amplitude [db]')
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

%% Simulink model with single controller

% Controller for mass 2
controller = 2;
switch controller
    case 1 % 23 Hz
        % Gain
        H_gain = 2;  
        % Integrator
        I_I = 5;
        H_I = (1+(I_I*2*pi/s)); 
        % Notch
        Nf1 = 58.1;     % Location of zeros [hz]
        Nf2 = 80;     % Location of poles [hz]
        Nbeta1 = 0.02;
        Nbeta2 = 5;
        H_notch = ((1/(2*pi*Nf1)^2)*s^2+(2*Nbeta1/(2*pi*Nf1))*s+1)/((1/(2*pi*Nf2)^2)*s^2+(2*Nbeta2/(2*pi*Nf2))*s+1); 
        % PD
        P = 1;
        D = 0.025;
        H_pd = P+D*s;  
        % Lead/lag controller
        Lead_zero = 6.667;
        Lead_pole = 90;
        H_Lead = ((s/(2*pi*Lead_zero))+1)/((s/(2*pi*Lead_pole))+1);   
        % lowpass
        f_LP = 550;
        H_LP = 1/((1/(2*pi*f_LP))*s+1);     
        % combine blocks
        C2 = H_gain*H_notch*H_Lead*H_pd*H_I*H_LP;
    case 2 % 45-26 Hz bw
        % Gain
        H_gain = 6;       
        % Integrator
        I_I = 13;
        H_I = (1+(I_I*2*pi/s));   
        % Notch
        Nf1 = 58;     % Location of zeros [hz]
        Nf2 = 150;     % Location of poles [hz]
        Nbeta1 = 0.02;
        Nbeta2 = 0.5;
        H_notch = ((1/(2*pi*Nf1)^2)*s^2+(2*Nbeta1/(2*pi*Nf1))*s+1)/((1/(2*pi*Nf2)^2)*s^2+(2*Nbeta2/(2*pi*Nf2))*s+1);  
        % PD
        P = 0.2;
        D = 0.01;
        H_pd = P+D*s;       
        % Lead/lag 
        Lead_zero = 80;
        Lead_pole = 130;
        H_Lead = ((s/(2*pi*Lead_zero))+1)/((s/(2*pi*Lead_pole))+1);  
        % Lowpass
        f_LP = 300;
        H_LP = 1/((1/(2*pi*f_LP))*s+1);        
        % combine blocks
        C2 = H_gain*H_notch*H_Lead*H_pd*H_I*H_LP; 
end        
        
C2_resp = squeeze(freqresp(C2,f*2*pi));

fdiff_2LP = 50;
Hdiff = 1/((1/(2*pi*fdiff_2LP))*s^2+(20/(2*pi*fdiff_2LP))*s+1);
bode(Hdiff)

% Notch filter for FF
NotchFilter = true
switch NotchFilter
    case (true)
        NFFf1 = 58;     % Location of zeros [hz]
        NFFf2 = 58;     % Location of poles [hz]
        NFFbeta1 = 0.02;
        NFFbeta2 = 1;
        CFF = ((1/(2*pi*NFFf1)^2)*s^2+(2*NFFbeta1/(2*pi*NFFf1))*s+1)/((1/(2*pi*NFFf2)^2)*s^2+(2*NFFbeta2/(2*pi*NFFf2))*s+1);
    case (false)
        CFF = tf([1],[1]);
end

if (true)
    figure(4)
    subplot(2,1,1)
    semilogx(f,db(C2_resp),f,db(H2_model),f,db(C2_resp.*H2_model));
    grid on;
    title('System from motor to mass 1')
    ylabel('Amplitude [db]')
    legend('C1','H1','C1*H1')
    subplot(2,1,2)
    semilogx(f,(angle(C2_resp)*180/pi),f,(angle(H2_model)*180/pi),f,(angle(C2_resp.*H2_model)*180/pi));
    grid on;
    ylabel('Angle [deg]')
    xlabel('Frequency [hz]')
end


%% FRF from simulation to check if everyting is correct

SimTime = 5; % [sec]

disp('Starting simulation')
sim('simulation_singlecontroller');
disp('Simulation finished')

nfft_sim = 5*(1/SampleTime/10); % Window length

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


