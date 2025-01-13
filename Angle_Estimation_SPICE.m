%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numElements = 10    ;                        % Number of array elements
numSnapshots = 500;                  % Number of snapshots
elementSpacing = 0.5;                 % Spacing between array elements in wavelengths                                  % SNR values in db
omg = [10 -5];
N = numElements^2;  N1=N^2;% Number of Samples
numSignals = 2;                          % Number of incoming signals
SNR = 10;                                     % SNR values in db
M = 512;      M1 = M+N;       Im = eye(N);    N2 = sqrt(N);                                      % Extended Number of Samples for high resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the steering vector and generating the incoming signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



rad_angles = deg2rad(omg);
eomg = linspace(-90,90,M);

steeringVector = exp(-1j*(0:numElements-1)'*pi*sin(rad_angles));


signals = (randn(numSignals, numSnapshots) + 1i * randn(numSignals, numSnapshots)) * sqrt(0.5); % Generating a Pseudo-Random Signal


% Generate noise
noise = (randn(numElements, numSnapshots) + 1i * randn(numElements, numSnapshots)) * sqrt(0.5);  % Generation of Pseudo-Random Noise
noisePower = 10^(-SNR/10); % Noise Power
noise = sqrt(noisePower) * noise;

% Received Signal Generation
R = steeringVector*signals + noise;
R_signal = (R*R')/numSnapshots;
singals = reshape(R_signal,[],1); % vectorized signal correlation matrix


%% Power Computation
steeringMatrix = exp(-1j*pi*(0:numElements-1)'*sin(deg2rad(eomg)));
DFT=[];
for i=1:M
    dd = steeringMatrix(:,i);  % Column vector of steering matrix
    ddh = dd*dd';                   % Computation of hermitian product
    ddhh = reshape(ddh,[],1);  % vectorization of hermitian product
    DFT = [DFT,ddhh];                           % Adding the vectorized product vectors of all columns
    clear dd;clear ddh;

end
DFT = [DFT Im];
IDFT = DFT';



%%
% Finding the Normalized Power
Signal_Power = (abs(fft(singals,M))).^2;  % Computation of FFT
normsignal_Power = Signal_Power/(N1);  % Dividing by norm of Steering Vector
S_Po=[normsignal_Power; ((abs(singals).^2))];




N_S=norm(singals);

w=[(N2/N_S)*ones(1,M) (1/N_S)*ones(1,(M1-M))];

s_old= S_Po;
for jj=1:1e6
    
    R = DFT*diag(s_old)*IDFT;
    R_inv=inv(R);
    u = R_inv*singals;
    c = abs(IDFT*u);

    rho = (w*(s_old.*c));

    s_new = (1/rho)*(s_old.*c).*(1./w');


    Er_En = (norm(s_old-s_new))/(norm(s_old));

    if(Er_En>1e-5)

        s_old = s_new;
    else
        break;
    end
end
s_new=real(s_new);
snew1 = s_new(1:M);



figure; plot(eomg,snew1);
title("Power Spectrum Using SPICE Optimization");
xlabel("Normalized Frequency");ylabel("Amplitude");
