function u = waveletOpt(Z,lambda2,options)
% WAVELETOPT: find optimal orthogonal wavelet basis. 
% 
%   U = WAVELETOPT(Z,LAMBDA2,OPTIONS) returns the optimal wavelet based on 
%   the data in Z.
% 
%   INPUTS
%       Z           data matrix (m-by-n, assuming m divisible by 2)
%       lambda2     pre-factor for variance penalty
%       options     options for fminunc
% 
%   OUTPUTS
%       u           wavelet (m-by-1)
%
%   NOTES
%       none
%
%
%   FLORYAN, DANIEL
%   February 29, 2020
%   Edited July 17, 2020

[m,n] = size(Z);
uhat = zeros(m,1);

% Set initial guess for gamma and theta
if mod(m,4)==0 % m divisible by 4
    gamma0 = 2*pi*rand(m/4,1); % random initial guess for gamma
    theta0 = 2*pi*rand(m/2 - 1,1); % random initial guess for theta
else % m not divisible by 4
    gamma0 = 2*pi*rand((m-2)/4 + 1,1); % random initial guess for gamma
    theta0 = 2*pi*rand(m/2 - 1,1); % random initial guess for theta
end

% Perform optimization
ZZT = Z*Z.'; % pre-compute this once
x = fminunc(@(x) costfun(x,Z,ZZT,lambda2),[gamma0;theta0],options);

% Extract gamma and theta, construct r and u
if mod(m,4)==0
    gamma = x(1:m/4);
    theta = x(m/4+1:end);
    r = [sqrt(2)*cos(gamma); 1; sqrt(2)*sin(flipud(gamma))];
else
    gamma = x(1:(m-2)/4+1);
    theta = x((m-2)/4+2:end);
    r = [sqrt(2)*cos(gamma); sqrt(2)*sin(flipud(gamma))];
end
uhat(1:m/2+1) = [r(1); r(2:m/2).*exp(1i*theta); r(m/2+1)];
uhat(m/2+2:m) = conj(flipud(uhat(2:m/2)));
u = ifft(uhat);


function [f,dfdx] = costfun(x,Z,ZZT,lambda2)
[m,n] = size(Z);
uhat = zeros(m,1);

% Extract gamma and theta, construct r and u
if mod(m,4)==0
    gamma = x(1:m/4);
    theta = x(m/4+1:end);
    r = [sqrt(2)*cos(gamma); 1; sqrt(2)*sin(flipud(gamma))];
else
    gamma = x(1:(m-2)/4+1);
    theta = x((m-2)/4+2:end);
    r = [sqrt(2)*cos(gamma); sqrt(2)*sin(flipud(gamma))];
end
uhat(1:m/2+1) = [r(1); r(2:m/2).*exp(1i*theta); r(m/2+1)];
uhat(m/2+2:m) = conj(flipud(uhat(2:m/2)));
u = real(ifft(uhat));

% Calculate cost
% Construct U differently based on size (to increase speed)
if m > 64
    U = zeros(m,m);
    for i=1:m
        U(:,i) = circshift(u,i-1);
    end
else
    U = real(ifft(diag(uhat))*fft(eye(m)));
end
f = -sum(sum((Z.'*U(:,1:2:m)).^2))/sum(sum(Z.^2));

% Calculate variance and add to cost function
xbar = sum(cos(2*pi*[0:m-1]'/m).*u.*u);
ybar = sum(sin(2*pi*[0:m-1]'/m).*u.*u);
rbar = sqrt(xbar*xbar + ybar*ybar);
uvar = 1 - rbar;
f = f + lambda2*uvar;

% Calculate gradient of variance wrt u (will be used later)
duvardu = lambda2*2/(uvar-1)*(xbar*cos(2*pi*[0:m-1]/m) + ybar*sin(2*pi*[0:m-1]/m)).*u';

% Calculate gradient of cost function wrt gamma and theta
dfdx = zeros(length(x),1);
if m > 1024
    temp = bsxfun(@times,fft(ZZT),conj(uhat));
    temp = ifft(temp);
    temp(2:2:end,:) = 0;
    temp = fft(temp);
    temp = ifft(temp.'); % this matrix will be used for all derivatives
else
    temp = U;
    temp(:,2:2:end) = 0;
    temp = ZZT*temp; % this order of multiplication is faster assuming n>m
    temp = fft(temp.');
    temp = ifft(temp.'); % this matrix will be used for all derivatives
end
Znorm2 = sum(sum(Z.^2));

if mod(m,4)==0
    % Derivative of cost wrt gamma_0
    dfdx(1) = -(-temp(1,1)*sqrt(2)*sin(gamma(1)) + ...
               temp(m/2 + 1,m/2 + 1)*sqrt(2)*cos(gamma(1)))*2/Znorm2 + ...
               duvardu*ifft([-sqrt(2)*sin(gamma(1)); zeros(m/2-1,1); sqrt(2)*cos(gamma(1)); zeros(m/2-1,1)]);
           
    % Derivative of cost wrt gamma_j for j=1,...,m/4 - 1
    for j=1:(m/4 - 1)
        dfdx(1+j) = -(-temp(j+1,j+1)*sqrt(2)*sin(gamma(j+1))*exp(1i*theta(j)) + ...
                     temp(m/2 - j+1,m/2 - j+1)*sqrt(2)*cos(gamma(j+1))*exp(1i*theta(m/2 - j)) + ...
                     temp(m/2 + j+1,m/2 + j+1)*sqrt(2)*cos(gamma(j+1))*exp(-1i*theta(m/2 - j)) + ...
                    -temp(m-j+1,m-j+1)*sqrt(2)*sin(gamma(j+1))*exp(-1i*theta(j)))*2/Znorm2 + ...
                     duvardu*ifft([zeros(j,1); -sqrt(2)*sin(gamma(j+1))*exp(1i*theta(j)); zeros(m/2-2*j-1,1); sqrt(2)*cos(gamma(j+1))*exp(1i*theta(m/2 - j)); zeros(2*j-1,1); sqrt(2)*cos(gamma(j+1))*exp(-1i*theta(m/2 - j)); zeros(m/2-2*j-1,1); -sqrt(2)*sin(gamma(j+1))*exp(-1i*theta(j)); zeros(j-1,1)]);
    end
    
    % Derivative of cost wrt theta_j for j=1,...,m/2 - 1
    for j=1:(m/2 - 1)
        dfdx(m/4 + j) = -(temp(j+1,j+1)*uhat(j+1)*1i + ...
                        -temp(m-j+1,m-j+1)*uhat(m-j+1)*1i)*2/Znorm2 + ...
                        duvardu*ifft([zeros(j,1); uhat(j+1)*1i; zeros(m-2*j-1,1); -uhat(m-j+1)*1i; zeros(j-1,1)]);
    end
else
    % Derivative of cost wrt gamma_0
    dfdx(1) = -(-temp(1,1)*sqrt(2)*sin(gamma(1)) + ...
               temp(m/2 + 1,m/2 + 1)*sqrt(2)*cos(gamma(1)))*2/Znorm2 + ...
               duvardu*ifft([-sqrt(2)*sin(gamma(1)); zeros(m/2-1,1); sqrt(2)*cos(gamma(1)); zeros(m/2-1,1)]);
        
    % Derivative of cost wrt gamma_j for j=1,...,(m-2)/4
    for j=1:(m-2)/4
        dfdx(1+j) = -(-temp(j+1,j+1)*sqrt(2)*sin(gamma(j+1))*exp(1i*theta(j)) + ...
                     temp(m/2 - j+1,m/2 - j+1)*sqrt(2)*cos(gamma(j+1))*exp(1i*theta(m/2 - j)) + ...
                     temp(m/2 + j+1,m/2 + j+1)*sqrt(2)*cos(gamma(j+1))*exp(-1i*theta(m/2 - j)) + ...
                    -temp(m-j+1,m-j+1)*sqrt(2)*sin(gamma(j+1))*exp(-1i*theta(j)))*2/Znorm2 + ...
                    duvardu*ifft([zeros(j,1); -sqrt(2)*sin(gamma(j+1))*exp(1i*theta(j)); zeros(m/2-2*j-1,1); sqrt(2)*cos(gamma(j+1))*exp(1i*theta(m/2 - j)); zeros(2*j-1,1); sqrt(2)*cos(gamma(j+1))*exp(-1i*theta(m/2 - j)); zeros(m/2-2*j-1,1); -sqrt(2)*sin(gamma(j+1))*exp(-1i*theta(j)); zeros(j-1,1)]);
    end
    
    % Derivative of cost wrt theta_j for j=1,...,m/2 - 1
    for j=1:(m/2 - 1)
        dfdx((m-2)/4 + 1 + j) = -(temp(j+1,j+1)*uhat(j+1)*1i + ...
                                -temp(m-j+1,m-j+1)*uhat(m-j+1)*1i)*2/Znorm2 + ...
                                duvardu*ifft([zeros(j,1); uhat(j+1)*1i; zeros(m-2*j-1,1); -uhat(m-j+1)*1i; zeros(j-1,1)]);
    end
end

dfdx = real(dfdx);
