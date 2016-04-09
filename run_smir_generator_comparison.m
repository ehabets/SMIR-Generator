%% Initialization
close all
clear
clc

%% Setup
procFs = 8000;                      % Sampling frequency (Hz)
c = 343;                            % Sound velocity (m/s)
nsample = 512;                      % Length of desired RIR
N_harm = 20;                        % Maximum order of harmonics to use in SHD
K = 2;                              % Oversampling factor

L = [4 6 8];                        % Room dimensions (x,y,z) in m
sphLocation = [2 3.2 4];            % Receiver location (x,y,z) in m
s = [2.37 4.05 4.4];                % Source location(s) (x,y,z) in m
beta = [1 0.7 0.7 0.5 0.2 1];       % Room reflection coefficients [\beta_x_1 \beta_x_2 \beta_y_1 \beta_y_2 \beta_z_1 \beta_z_2]
order = -1;                         % Reflection order (-1 is maximum reflection order)

sphRadius = 0.042;                  % Radius of the sphere (m)
sphType = 'open';                   % Type of sphere (open/rigid)

mic = [pi/4 pi/4; pi/2 pi/4];	    % Microphone positions (azimuth, inclination)

src_type = 'o';                     % Directional source type ('o','c','s','h','b')
[src_ang(1),src_ang(2)] = mycart2sph(sphLocation(1)-s(1),sphLocation(2)-s(2),sphLocation(3)-s(3)); % Towards the receiver

%% Run spherical array simulation
[h, H, beta_hat] = smir_generator(c, procFs, sphLocation, s, L, beta, sphType, sphRadius, mic, N_harm, nsample, K, order, 0, 0, src_type, src_ang);

%% Compare with RIR generator run separately for each microphone in the array
[mic_pos(:,1), mic_pos(:,2), mic_pos(:,3)] = mysph2cart(mic(:,1),mic(:,2),sphRadius); % Microphone positions relative to centre of array ("sphLocation")
mic_pos = mic_pos + repmat(sphLocation,size(mic,1),1);

h_rirgen = 4*pi*rir_generator(c, procFs, mic_pos, s, L, beta, nsample, 'omnidirectional', order, 3, [0 0], false);

H_rirgen = fft(h_rirgen, [], 2);

%% Plotting

mic_to_plot = 1;

figure;
subplot(211);
plot([0:nsample-1]/procFs, h_rirgen(mic_to_plot,1:nsample), 'g')
hold all;
plot([0:nsample-1]/procFs,h(mic_to_plot,1:nsample), 'r')
xlim([0 (nsample-1)/procFs]);
title(['Room impulse response at microphone ', num2str(mic_to_plot)]);
xlabel('Time (s)');
ylabel('Amplitude');
legend('RIR generator', 'SMIR generator');

subplot(212);
plot((0:1/nsample:1/2)*procFs,mag2db(abs(H_rirgen(mic_to_plot,1:nsample/2+1))), 'g');
hold all;
plot((0:1/(K*nsample):1/2)*procFs,mag2db(abs(H(mic_to_plot,1:K*nsample/2+1))), 'r');
title(['Room transfer function magnitude at microphone ', num2str(mic_to_plot)]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
legend('RIR generator', 'SMIR generator');