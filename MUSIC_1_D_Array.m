%%
 % Copyright (c) 2025, Sanjeeva Reddy S
 % All rights reserved.
 
 %This source code is licensed under the MIT license found in the
 % LICENSE file in the root directory of this source tree.
 
 % Unauthorized copying of this file, via any medium, is strictly prohibited
 % unless explicit permission is granted by the copyright owner.
 
 % Description:
 % This file contains utility functions for processing sparse arrays.
 
 % Author: Sanjeeva Reddy S
 % EMail: sanjeevareddy.s414@gmail.com
 % Created on: January 5, 2025



%% Simulation of 1-Dimensional Array for Estimation of Angle of Arrival Using MUSIC Algorithm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close ;tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numElements = 8;                        % Number of array elements
numSignals = 3;                           % Number of incoming signals
numSnapshots = 100;                  % Number of snapshots
elementSpacing = 0.5;                 % Spacing between array elements in wavelengths
SNR = 0;                                    % SNR values in db
angles = [50 70 25];                     % Angles of arrival in degrees


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


R = steeringVector*signals + noise;
R_signal = (R*R')/numSnapshots;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of Angle of Arrival using MUSIC Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eigen Value Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Perform Eigenvalue Decomposition
[eigenVectors, eigenValues] = eig(R_signal);  % Eigen value computation for covariance matrix
[~, idx] = sort(diag(eigenValues), 'descend');  
eigenVectors = eigenVectors(:, idx);  % Signal Subspace

% Estimate noise subspace
noiseSubspace = eigenVectors(:, numSignals+1:end); % Last N-numSignals eigenvectors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spectrum Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = -90:1:90;

musicSpectrum = zeros(length(theta), 1);  % Music Spectrum Matrix

% Calculate the MUSIC spectrum
a = exp(-1j*pi*(0:numElements-1)'*sin(deg2rad(theta)));
for k = 1:length(theta)
    musicSpectrum(k) = 1 / norm(noiseSubspace' * a(:,k), 2)^2; % Formula for Spectrum Calculation
end


% Normalize the MUSIC spectrum
musicSpectrum = 10 * log10(musicSpectrum / max(musicSpectrum)); % Normalizing the spectrum values


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the Spectrum and Angle Estimation from Plot


% Plot the MUSIC spectrum
figure;
plot(theta, musicSpectrum);
xlabel('Angle (degrees)');
ylabel('MUSIC Spectrum (dB)');
title('MUSIC Spectrum');
grid on;

% Find the peaks in the MUSIC spectrum
[~, peakIndices] = findpeaks(musicSpectrum, 'SortStr', 'descend', 'NPeaks', numSignals);
estimatedAngles = theta(peakIndices);

% Display the estimated angles
disp('Estimated DOAs:');
disp(estimatedAngles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
