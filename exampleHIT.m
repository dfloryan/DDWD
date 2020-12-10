% Script to find and plot optimal wavelets for HIT data. 
% Roughly reproduces figure 6 in the paper.
% WARNING: will take a long time to run. 

% Load HIT dataset
load('data/isotropic4096long.mat')
Z = ulong;
[m,n] = size(Z);

% Set some parameters
levels = 12; % = log2(m)
lambda2 = 10; % variance penalty

% Set options for gradient descent
options = optimoptions(@fminunc,'Algorithm','quasi-newton',...
    'CheckGradients',false,...
    'Display','iter-detailed',...
    'HessUpdate','bfgs',...
    'MaxIterations',5000,...
    'OptimalityTolerance',1e-6,...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance',1e-6,...
    'FunctionTolerance',1e-6);

% Pre-allocate some matrices
W = Z; % stores wavelet coefficients
u = zeros(m,levels) + NaN; % this will store all the filters

% Find the optimal filters, compute resulting wavelet coefficients
for i=1:levels
    u(1:m/2^(i-1),i) = waveletOpt(W(m-m/2^(i-1)+1:m,:),lambda2,options);
    W(m-m/2^(i-1)+1:m,:) = dwtos(W(m-m/2^(i-1)+1:m,:),u(1:m/2^(i-1),i),u2v(u(1:m/2^(i-1),i)));
    lambda2 = lambda2/4; % update variance penalty
end

% Plot the wavelets we computed, offset from each other by 0.25, as well as
% their power spectra
cols = copper(levels + 1);
for j=1:levels+1
    ind = 2^levels - 2^(j-1) + 1;
    Zr = zeros(m,1);
    Zr(ind) = 1;
    for i=levels:-1:1
        Zr(m-m/2^(i-1)+1:m,:) = idwtos(Zr(m-m/2^(i-1)+1:m,:),u(1:m/2^(i-1),i),u2v(u(1:m/2^(i-1),i)));
    end
    
    figure(1)
    hold on
    plot([0:m-1],(j-1)*0.25 + Zr,'-','color',cols(j,:),'linewidth',2)
    
    figure(2)
    Zpow = abs(fft(Zr)).^2;
    Zpow = Zpow(1:m/2 + 1);
    Zpow(2:m/2) = 2*Zpow(2:m/2);
    hold on
    plot(0:m/2,Zpow/m,'-','color',cols(j,:),'linewidth',2)
end
figure(1)
xlabel('j')
figure(2)
xlabel('k')
