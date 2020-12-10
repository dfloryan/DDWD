function z = idwtos(w,u,v)
% IDWTOS: inverse one-stage discrete wavelet transform.
% 
%   Z = IDWTOS(W,U,V) returns the signal generated by one-stage wavelet 
%   coefficients w with filters u and v.
% 
%   INPUTS
%       w           wavelet coefficients (N-by-m, assuming N = 2*M, m
%                   signals)
%       u           low-pass filter (N-by-1)
%       v           high-pass filter (N-by-1)
% 
%   OUTPUTS
%       z           signal (N-by-m)
%
%   NOTES
%       conventions as in Frazier book
%
%
%   FLORYAN, DANIEL
%   May 1, 2020
%   Edited May 13, 2020

[N,m] = size(w);

% Need to handle N==2 separately because of how upsample function works
if N==2
    sig1 = [w(1:N/2,:); zeros(1,m)];
    sig2 = [w(N/2+1:N,:); zeros(1,m)];
else
    sig1 = upsample(w(1:N/2,:),2);
    sig2 = upsample(w(N/2+1:N,:),2);
end
z = ifft(bsxfun(@times,fft(sig1),fft(v))) + ifft(bsxfun(@times,fft(sig2),fft(u)));
