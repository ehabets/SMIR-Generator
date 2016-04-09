%% Initialization
close all
clear
clc

%% Setup
procFs = 8000;                      % Sampling frequency (Hz)
c = 343;                            % Sound velocity (m/s)
nsample = 2*1024;                   % Length of desired RIR
N_harm = 30;                        % Maximum order of harmonics to use in SHD
K = 1;                              % Oversampling factor

L = [5 6 4];                        % Room dimensions (x,y,z) in m
sphLocation = [1.6 4.05 1.7];       % Receiver location (x,y,z) in m
s = [3.37 4.05 1.7];                % Source location(s) (x,y,z) in m

HP = 1;                             % Optional high pass filter (0/1)
src_type = 'o';                     % Directional source type ('o','c','s','h','b')
[src_ang(1),src_ang(2)] = mycart2sph(sphLocation(1)-s(1),sphLocation(2)-s(2),sphLocation(3)-s(3)); % Towards the receiver

%% Example 1
order = 6;                          % Reflection order (-1 is maximum reflection order)
refl_coeff_ang_dep = 0;             % Real reflection coeff(0) or angle dependent reflection coeff(1)
beta = 0.3;                         % Reverbration time T_60 (s)
% beta = 0.2*ones(1,6);             % Room reflection coefficients [\beta_x_1 \beta_x_2 \beta_y_1 \beta_y_2 \beta_z_1 \beta_z_2]

sphRadius = 0.042;                  % Radius of the spherical microphone array (m)
sphType = 'rigid';                  % Type of sphere (open/rigid)
mic = [pi/4 pi; pi/2 pi];		    % Microphone positions (azimuth, elevation)

[h1, H1] = smir_generator(c, procFs, sphLocation, s, L, beta, sphType, sphRadius, mic, N_harm, nsample, K, order, refl_coeff_ang_dep, HP, src_type, src_ang);

%% Example 2
order = 6;                          % Reflection order (-1 is maximum reflection order)
refl_coeff_ang_dep = 1;             % Real reflection coeff(0) or angle dependent reflection coeff(1)
sigma = 1.5*10^4*ones(1,6);         % "Effective" flow resistivity; typical values between 10^3 and 10^9

sphRadius = 0.042;                  % Radius of the spherical microphone array (m)
sphType = 'rigid';                  % Type of sphere (open/rigid)
mic = [pi/4 pi; pi/2 pi];		    % Microphone positions (azimuth, elevation)

[h2, H2] = smir_generator(c, procFs, sphLocation, s, L, sigma, sphType, sphRadius, mic, N_harm, nsample, K, order, refl_coeff_ang_dep, HP, src_type, src_ang);

%% Example 3 (single microphone)

order = 6;                          % Reflection order (-1 is maximum reflection order)
refl_coeff_ang_dep = 0;             % Real reflection coeff(0) or angle dependent reflection coeff(1)
beta = 0.3;                         % Reverbration time T_60 (s)
% beta = 0.2*ones(1,6);             % Room reflection coefficients [\beta_x_1 \beta_x_2 \beta_y_1 \beta_y_2 \beta_z_1 \beta_z_2]

sphRadius = 0;                      % Single microphone at the center of the sphere with radius zero
sphType = 'open';                   % Type of sphere (open/rigid)
mic = [0 0];			            % Microphone positions (azimuth, elevation)

[h3, H3] = smir_generator(c, procFs, sphLocation, s, L, beta, sphType, sphRadius, mic, N_harm, nsample, K, order, refl_coeff_ang_dep, HP, src_type, src_ang);

%% Plotting

mic_to_plot = 1;

figure(1);
ax1(1)=subplot(211);
plot([0:nsample-1]/procFs,h1(mic_to_plot,1:nsample), 'r')
xlim([0 (nsample-1)/procFs]);
title(['Room impulse response at microphone ', num2str(mic_to_plot),' (real refl coeff)']);
xlabel('Time (s)');
ylabel('Amplitude');

ax1(2)=subplot(212);
plot([0:nsample-1]/procFs,h2(mic_to_plot,1:nsample), 'r')
xlim([0 (nsample-1)/procFs]);
title(['Room impulse response at microphone ', num2str(mic_to_plot), ' (angle dependent refl coeff)']);        % open sphere
xlabel('Time (s)');
ylabel('Amplitude');
linkaxes(ax1);