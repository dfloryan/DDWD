% Script to find and plot optimal wavelets for Gaussian white noise data. 
% Roughly reproduces figure 3 in the paper. 

% Create dataset of Gaussian white noise
Z = randn(2^5,1e5);
[m,n] = size(Z);

% Set some parameters
levels = 5; % = log2(m)
lambda2 = 0.02; % variance penalty

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

% Plot the wavelets we computed, offset from each other by 1
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
    plot([0:m-1],(j-1)*1 + Zr,'-','color',cols(j,:),'linewidth',2)
end
xlabel('j')
