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





%% Simulation of 1-Dimensional Array for Estimation of Angle of Arrival Using ROOT - MUSIC Algorithm for Coherent Sources

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;close all;tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numElements = 8;% Number of antenna elements
elementSpacing = 0.5;% Half wavelength spacing
numSignals = 3; % Number of signal sources
theta = [30 45 60]; % Directions of arrival
numSnapshots = 100;% Number of samples
signalFreq = 1e9;                         % Frequency of incoming signals (Hz)
c = 3e8;                                        % Speed of light (m/s)
wavelength = c / signalFreq;        % Wavelength of incoming signals (m)
SNR = 20;                                  % Signal to Noise Ratio (dB)
q = numElements-numSignals;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the steering vectors and Singal matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Construct steering matrix A for all sources
SteeringVector = [];
for n = 1:numSignals
    a = exp(-1i * 2 * pi * elementSpacing * sin(theta(n) * pi / 180) * (0:numElements-1)');
    SteeringVector = [SteeringVector, a];
end


% Generate signal sources
t = 1:numSnapshots;  % Time samples

% Coherent signals: same signal, but two of them are coherent (reflections)
signal1 = exp(1j * 2 * pi * 0.01 * t);  % First coherent signal
signal2 = signal1;                      % Second coherent signal (coherent with signal1)
signal3 = exp(1j * 2 * pi * 0.02 * t);  % Third non-coherent signal

signals = [signal1;signal2;signal3];


% Generating noise signal

noise_variance = 10^(-SNR/10);
noise = sqrt(noise_variance/2) * (randn(numElements, numSnapshots) + 1i * randn(numElements, numSnapshots));

R = SteeringVector*signals+ noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Applying Spatial Smoothening to Received Signal Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_signal = zeros(q,q);

for i=0:numElements-q
    z = (R((i+1:q+i),:));
    r = z*z'/(4*numSnapshots);
    R_signal = R_signal+r;

end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROOT MUSIC Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eigen Value Decomposition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Eigenvalue decomposition and sorting
[eigenVector, eigenValues] = eig(R_signal);
signalSubspace = zeros(1, q);


for k = 1:q
    signalSubspace(k) = real(eigenValues(k, k));
end
[signalSubspace, Ind] = sort(signalSubspace, 'descend');
signalSubspace = diag(signalSubspace);

noiseSubspace = zeros(q, q);
for k = 1:q
    noiseSubspace(:, k) = eigenVector(:, Ind(k));
end

% Constructing the noise subspace
En = [];
for d = numSignals+1:q
    En = [En, noiseSubspace(:, d)];
end

% Finding zeros of polynomial, choose the closest to the unit circle
C = En * En';
for k = (q-1):-1:-(q-1)
    P(q-k) = sum(diag(C, k));
end

% Computation of Roots of the Polynomial
rts = roots(P);
figure(2);zplane(rts, 1);title("Roots of the Ploynomial");


% Choose the roots closest to the unit circle
dist_from_unit_circle =abs(abs(rts - 1)); % Computing the roots nearer to the Unit Circle
theta = -180 * asin(angle(rts) / 2 / pi / elementSpacing) / pi;
disp("Estimated Angle of Arrivals: "+theta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
