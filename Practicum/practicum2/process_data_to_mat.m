%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%   4CM00 - Control Engineering                 %
% Practicum - Split measured data into a struct %
%                                               %
%   Author: Marcel van Wensveen                 %
%       &   Job Meijer                          %
%   Date:   23-10-2019                          %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all; close all; clc;
%% 

path = 'C:\Users\marce\Google Drive\S&C master\Job en marcel\4CM00 - control engineering\Practicum\Data\';
filename = 'extra';
load('20181023_1200.mat')
tempData = ans';

% 1. Time
% 2. enc0 [y1]
% 3. enc1 [y2]
% 4. control output (dac input [u]
% 5. error (control input) [e]
% 6. disturbance input (inserted at control input) [d]
% 7. input/reference [r]

Data.Time = tempData(:,1);
Data.y1 = tempData(:,2);
Data.y2 = tempData(:,3);
Data.u = tempData(:,4);
Data.e = tempData(:,5);
Data.d = tempData(:,6);
Data.r = tempData(:,7);

figure(3)
plot(Data.Time, Data.d);
% plot(Data.Time, Data.r,Data.Time, Data.y2);

save([path,filename,'.mat'], 'Data');