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




%% Simulation of 1-Dimensional Array for Estimation of Angle of Arrival Using ROOT - MUSIC Algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;close;tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numElements = 8;% Number of antenna elements
elementSpacing = 0.5;% Half wavelength spacing
numSignals = 3; % Number of signal sources
angles = [-20 45 20 ]; % Directions of arrival
numSnapshots = 200;% Number of samples
SNR = 20;                                  % Signal to Noise Ratio (dB)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the steering vectors and Singal matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steeringVector = exp(-1j*(0:numElements-1)'*pi*sin(deg2rad(angles)));


% Generating source signal
signal = randn(numSignals, numSnapshots);


% Generating noise signal

noise_variance = 10^(-SNR/10);
noise = sqrt(noise_variance/2) * (randn(numElements, numSnapshots) + 1i * randn(numElements, numSnapshots));


% Received signal matrix
R = steeringVector * signal + noise;

R_signal = (R*R')/numSnapshots;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROOT MUSIC Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Perform Eigenvalue Decomposition
[eigenVectors, eigenValues] = eig(R_signal);  % Eigen value computation for covariance matrix
[~, idx] = sort(diag(eigenValues));  
noiseSubspace = eigenVectors(:, idx(1:numElements-numSignals));  % Noise Subspace

% Finding zeros of polynomial, 
C = noiseSubspace * noiseSubspace';
for k = (-numElements+1):(numElements-1)
    P(numElements+k) = sum(diag(C, k));
    % P(numElements+k) = P(numElements+k)*exp(-1i*pi*numElements+k*sin(angle(P(numElements+k))))
end

% Computation of Roots of the Polynomial
rts = roots(P);

figure(1);zplane(rts);title("Roots of the Ploynomial");


% Choose the roots closest to the unit circle
dist_from_unit_circle =abs(abs(rts)-1); % Computing the roots nearer to the Unit Circle
figure(2);zplane(rts(dist_from_unit_circle<0.09));title("Roots of the Ploynomial nearer to unit circle");

angles = -rad2deg(asin(angle(rts(dist_from_unit_circle<0.09)) / pi));
disp("Estimated Angle of Arrivals: "+angles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
