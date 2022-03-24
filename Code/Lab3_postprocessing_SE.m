% This script plots ... 
% ----------------------------------------------------------------------
% Usage: 
% 
%
%  input:    ---------
%   ...
%
%  output:   ---------
%    several plots:
%     ...
%
% Author :  Stephanie Evers, stephaev@stud.ntnu.no
%           Norwegian University of Science and Technology
%           Trondheim
%
% Date   :  23 March 2022

clc
clear
close all

%% measurement data in form of peak to peak amplitudes
% Measuring the Frequency at Three Distinc Angles
phi_FR = [0 45 90];
freq_vec_FR = [8000:1000:16000]; % = [8 9 10 11 12 13 14 15 16];  % in kHz

peaks_0 = [2.09 2.12 2.77 3.12 3.19 ...
    3.05 2.79 2.295 1.455];         % peak to peak average at 0 degrees
peaks_45 = [0.855 1.085 1.33 1.75 2.165 ...
    2.385 2.35 1.985 1.42];         % peak to peak average at 45 degrees
peaks_90 = [1.455 1.76 2.17 2.485 2.59 ...
    2.61 2.35 2.035 1.625];         % peak to peak average at 90 degrees

%% measurement data in form of peak to peak amplitudes
% Measuring the Beam Pattern and Directivity at Three Distinct Frequencies
freq_vec_BP = [8000 12000 16000];
phi_BP = [-90:15:90];              % in degrees

rho8k  = [1.655 1.59 1.41 1.37 1.385 1.45 1.52 1.47 ...
    1.415 1.415 1.47 1.57 1.555];       % peak to peak average at 8 kHz
rho12k = [2.66 2.635 2.385 2.07 2.02 2.89 3.115 2.69 ...
    2.045 2.34 2.585 2.52 2.27];        % peak to peak average at 12 kHz
rho16k = [1.46 1.385 1.285 0.83 0.775 1.89 2.15 1.33 ...
    0.655 1.06 1.385 1.535 1.495];      % peak to peak average at 16 kHz

%% normalizing it and converting it into dB
max_tot = max([peaks_0 peaks_45 peaks_90 rho8k rho12k rho16k]);

% Frequency
peaks0_norm = db(peaks_0./max_tot);
peaks45_norm = db(peaks_45./max_tot);
peaks90_norm = db(peaks_90./max_tot);

% Beam Pattern
rho8k_norm  =  db(rho8k./max_tot);
rho12k_norm = db(rho12k./max_tot);
rho16k_norm = db(rho16k./max_tot); 

%% results for the frequency response at three angels
% three curves of different colors
figure(1)
semilogx(freq_vec_FR, peaks0_norm, ...
    freq_vec_FR, peaks45_norm, ...
    freq_vec_FR, peaks90_norm)
title('Frequency Response for three different angles')
xlabel('Freq [Hz]')
ylabel('Amplitude [dB]')
% xlim([f_l f_u]);
legend('0°', '45°', '90°', 'Location', 'southeast')
grid on

%% results of beam pattern = directivity (polar plots) for the three 
% frequencies, three curves of different colors
lim = -10;
NN = 2;

figure(2);
sth1 = polardb(phi_BP,rho12k_norm,lim,NN,'-*r')
hold on
sth2 = polardb(phi_BP,rho16k_norm,lim,NN,'-*g')
sth3 = polardb(phi_BP,rho8k_norm,lim,NN,'-*b')
hold off
legend('','','','','','','f = 12kHz','f = 16kHz','f = 8kHz','Location','best')
title('Beam Pattern')

%% computing Directivity D and Directivity Index DI 
% individually for left and right quadrant

% splitting into the two quadrants
phi_BP_I  = [0:15:90];              % in degrees
phi_BP_II = [-90:15:0];              % in degrees

rho8k_I  = rho8k_norm(:,7:13);       % peak to peak average at 8 kHz
rho12k_I = rho12k_norm(:,7:13);        % peak to peak average at 12 kHz
rho16k_I = rho16k_norm(:,7:13);

rho8k_II  = rho8k_norm(:,1:7);       % peak to peak average at 8 kHz
rho12k_II = rho12k_norm(:,1:7);        % peak to peak average at 12 kHz
rho16k_II = rho16k_norm(:,1:7);

% IS THIS CALCULATION CORRECT??????????????????????????????????????
% -> or do they have to be switched ? think logical how the measurement was
% proceeded

%% right quadrant = I
% D 
D_8k_I  = 2 ./ (sum((rho8k_I.^2) .* sind(phi_BP_I)));
D_12k_I = 2 ./ (sum((rho12k_I.^2) .* sind(phi_BP_I)));
D_16k_I = 2 ./ (sum((rho16k_I.^2) .* sind(phi_BP_I)));

% DI
DI_8k = 10*log10(D_8k_I);
DI_12k = 10*log10(D_12k_I);
DI_16k = 10*log10(D_16k_I);

%% left quadrant = II
% D
D_8khz_II  = 2 ./ (sum((rho8k_II.^2) .* sind(phi_BP_II)));
D_12khz_II = 2 ./ (sum((rho12k_II.^2) .* sind(phi_BP_II)));
D_16khz_II = 2 ./ (sum((rho16k_II.^2) .* sind(phi_BP_II)));

% DI
DI_8kh_II   = 10*log10(D_8khz_II);
DI_12khz_II = 10*log10(D_12khz_II);
DI_16khz_II = 10*log10(D_16khz_II);

% end of file