%% Simulation of 1-Dimensional Array for Estimation of Angle of Arrival Using MUSIC Algorithm for Coherent Sources


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close;tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numElements = 8;                        % Number of array elements
numSignals = 3;                           % Number of incoming signals
numSnapshots = 100;                  % Number of snapshots
elementSpacing = 1;                 % Spacing between array elements in wavelengths
signalFreq = 1e9;                         % Frequency of incoming signals (Hz)
c = 3e8;                                        % Speed of light (m/s)
wavelength = c / signalFreq;        % Wavelength of incoming signals (m)
SNR = 10;                                    % SNR values in db
angles = [30 45 60];                     % Angles of arrival in degrees
q = numElements-numSignals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the steering vectors matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steeringVector = exp(-1j*(0:numElements-1)'*pi*sin(deg2rad(angles)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the Incoming Signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate signal sources
t = 1:numSnapshots;  % Time samples

% Coherent signals: same signal, but two of them are coherent (reflections)
signal1 = exp(1j * 2 * pi * 0.01 * t);  % First coherent signal
signal2 = signal1;                      % Second coherent signal (coherent with signal1)
signal3 = exp(1j * 2 * pi * 0.02 * t);  % Third non-coherent signal

signals = [signal1;signal2;signal3];

% Generate noise
noise = (randn(numElements, numSnapshots) + 1i * randn(numElements, numSnapshots)) * sqrt(0.5);  % Generation of Pseudo-Random Noise
noisePower = 10^(-SNR/10); % Noise Power 
noise = sqrt(noisePower) * noise; 

% Received Signal Generation
R = steeringVector*signals+ noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Applying Spatial Smoothening for the Received Signal Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_signal = zeros(q,q);

for i=0:numElements-q
    z = (R((i+1:q+i),:));
    r = z*z'/(4*numSnapshots);
    R_signal = R_signal+r;

end
   


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

theta = (-90:1:90);
numAngles = length(theta);
musicSpectrum = zeros(numAngles, 1);  % Music Spectrum Matrix

% Calculate the MUSIC spectrum
a = exp(-1i*(0:q-1)'*pi*sin(deg2rad(theta)));
for k = 1:numAngles
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


