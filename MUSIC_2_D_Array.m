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



%% Simulation of 2-Dimensional Array for Estimation of Angle of Arrival Using MUSIC Algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the Cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;clear;close;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


frequency = 1e9;                   % Carrier frequency (1 GHz)
c = 3e8;                           % Speed of light (m/s)
wavelength = c / frequency;        % Wavelength
M = 4;                             % Number of antenna elements in the x-dimension
N = 4;                             % Number of antenna elements in the y-dimension
elementSpacing = wavelength / 2;   % Distance between elements (half-wavelength spacing)
numSignals = 3;                    % Number of sources
azimuths = [-50, 45,60];              % True azimuth angles of arrival (degrees)
elevations = [30, -45, 60];             % True elevation angles of arrival (degrees)
numSnapshots = 100;                % Number of snapshots



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the steering matrix for a URA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


steeringVector = zeros(M * N, numSignals);
for k = 1:numSignals
    theta_rad = deg2rad(azimuths(k));  % Convert azimuth to radians
    phi_rad = deg2rad(elevations(k));  % Convert elevation to radians
    
    for m1 = 1:M
        for m2 = 1:N
            m = (m1 - 1) * N + m2;  % Linear index for URA element (x, y)
            % Steering vector calculation (corrected)
            steeringVector(m, k) = exp(-1j * 2 * pi *((m1 - 1) * elementSpacing * sin(phi_rad) * cos(theta_rad) / wavelength + (m2 - 1) * elementSpacing * sin(phi_rad) * sin(theta_rad) / wavelength));
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Received Signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


S = randn(numSignals, numSnapshots) + 1j * randn(numSignals, numSnapshots);  % Source signals
X = steeringVector * S;  % Received signal at the array

% Add noise
noise = (randn(M * N, numSnapshots) + 1j * randn(M * N, numSnapshots)) * sqrt(0.5);
X = X + noise;

% Compute the covariance matrix
Rxx = (1 / numSnapshots) * (X * X');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Eigen Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[eigenVectors, eigenValues] = eig(Rxx);

% Sort the eigenvalues in descending order
[eigenValues_sorted, idx] = sort(diag(eigenValues), 'descend');
eigenVectors_sorted = eigenVectors(:, idx);

% Separate noise and signal subspaces
En = eigenVectors_sorted(:, numSignals + 1:end);  % Noise subspace


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MUSIC Spectrum Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define the search grid for azimuth and elevation
theta_scan = deg2rad(-90:1:90);  % Azimuth angle scan (in radians)
phi_scan = deg2rad(-90:1:90);    % Elevation angle scan (in radians)
P_music = zeros(length(phi_scan), length(theta_scan));

% MUSIC Spectrum calculation
for i = 1:length(phi_scan)
    for j = 1:length(theta_scan)
        theta_rad = theta_scan(j);
        phi_rad = phi_scan(i);
        a_theta_phi = zeros(M * N, 1);  % Steering vector for each angle pair
        
        for m1 = 1:M
            for m2 = 1:N
                m = (m1 - 1) * N + m2;  % Linear index for URA element (x, y)
                a_theta_phi(m) = exp(-1j * 2 * pi * ((m1 - 1) * elementSpacing * sin(phi_rad) * cos(theta_rad) / wavelength + (m2 - 1) * elementSpacing * sin(phi_rad) * sin(theta_rad) / wavelength));
            end
        end
        
        % MUSIC pseudo-spectrum
        P_music(i, j) = 1 / (a_theta_phi' * (En * En') * a_theta_phi);
    end
end

% Convert to dB scale
P_music_dB = 10 * log10(abs(P_music));

% Plot the MUSIC spectrum
figure;
mesh(rad2deg(theta_scan), rad2deg(phi_scan), P_music_dB);
xlabel('Azimuth (degrees)');
ylabel('Elevation (degrees)');
zlabel('Spectrum Magnitude (dB)');
title('2D MUSIC Spectrum');

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
