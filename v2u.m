function u = v2u(v)
% V2U: create low-pass filter from high-pass filter.
% 
%   U = V2U(V) returns the low-pass wavelet filter corresponding to the
%   given high-pass filter. 
% 
%   INPUTS
%       v           high-pass filter (N-by-1)
% 
%   OUTPUTS
%       u           low-pass filter (N-by-1)
%
%   NOTES
%       conventions as in Frazier book
%
%
%   FLORYAN, DANIEL
%   May 5, 2020
%   Edited May 5, 2020

N = length(v);
u = zeros(N,1);
for i=1:N
    u(mod(3-i-1,N)+1) = (-1)^i*conj(v(i));
end
