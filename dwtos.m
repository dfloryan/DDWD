function w = dwtos(z,u,v)
% DWTOS: one-stage discrete wavelet transform.
% 
%   W = DWTOS(Z,U,V) returns the one-stage wavelet coefficients of the 
%   signal z with filters u and v.
% 
%   INPUTS
%       z           input signal (N-by-m, assuming N = 2*M, m signals)
%       u           low-pass filter (N-by-1)
%       v           high-pass filter (N-by-1)
% 
%   OUTPUTS
%       w           wavelet coefficients (N-by-m)
%
%   NOTES
%       conventions as in Frazier book
%
%
%   FLORYAN, DANIEL
%   May 1, 2020
%   Edited May 5, 2020

% Compute conjugate reflections of filters
utilde = conj(circshift(flipud(u),1));
vtilde = conj(circshift(flipud(v),1));

% Convolve and downsample
x = downsample(ifft(bsxfun(@times,fft(z),fft(vtilde))),2);
y = downsample(ifft(bsxfun(@times,fft(z),fft(utilde))),2);

w = [x;y];
