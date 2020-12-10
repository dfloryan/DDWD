function v = u2v(u)
% U2V: create high-pass filter from low-pass filter.
% 
%   V = U2V(U) returns the high-pass wavelet filter corresponding to the
%   given low-pass filter. 
% 
%   INPUTS
%       u           low-pass filter (N-by-1)
% 
%   OUTPUTS
%       v           high-pass filter (N-by-1)
%
%   NOTES
%       conventions as in Frazier book
%
%
%   FLORYAN, DANIEL
%   May 5, 2020
%   Edited May 5, 2020

N = length(u);
v = zeros(N,1);
for i=1:N
    v(i) = (-1)^i*conj(u(mod(3-i-1,N)+1));
end
