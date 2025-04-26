%%  Simulation of 1-Dimensional Array for Estimation of Angle of Arrival Using Conventional Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close ;tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numElements = 32;                        % Number of array elements
numSignals = 3;                           % Number of incoming signals
numSnapshots = 300;                  % Number of snapshots
elementSpacing = 0.5;                 % Spacing between array elements in wavelengths
SNR = 0;                                    % SNR values in db
angles = [50 -70 25];                     % Angles of arrival in degrees


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the steering vectors matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steeringVector = exp(-1j*(0:numElements-1)'*pi*sin(deg2rad(angles)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the Incoming Signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signals = (randn(numSignals, numSnapshots) + 1i * randn(numSignals, numSnapshots)) * sqrt(0.5); % Generating a Pseudo-Random Signal


% Generate noise
noise = (randn(numElements, numSnapshots) + 1i * randn(numElements, numSnapshots)) * sqrt(0.5);  % Generation of Pseudo-Random Noise
noisePower = 10^(-SNR/10); % Noise Power 
noise = sqrt(noisePower) * noise; 


R1 = steeringVector*signals + noise;
R = (R1*R1')/numSnapshots;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOA Estimation using Delay-and-Sum Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_scan = -90:0.5:90;   % Angle range to scan [degrees]
P = zeros(size(theta_scan));  % Power for each angle
a = exp(-1j*pi*(0:numElements-1)'*sin(deg2rad(theta_scan)));
% Beamforming (Delay-and-Sum)
for i = 1:length(theta_scan)
    P(i) = abs(a(:,i)' * R * R' * a(:,i));  % Beamformer output power
end

% Normalize and plot the results
P = P / max(P);
figure;
plot(theta_scan, 10*log10(P));
xlabel('Angle (degrees)');
ylabel('Power (dB)');
title('DOA Estimation using Delay-and-Sum Beamforming');
