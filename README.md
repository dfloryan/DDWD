# DDWD
Data-driven wavelets

README for Data-Driven Wavelet Decomposition (DDWD)

This distribution contains the code needed to create data-driven wavelets, 
as described in "Discovering multiscale and self-similar structure with 
data-driven wavelets," by D. Floryan and M. D. Graham, PNAS, 2020. 

This distribution contains five primary MATLAB functions: 
(1) waveletOpt.m: calculates data-driven wavelets
(2) dwtos.m: one-stage discrete wavelet transform
(3) idwtos.m: inverse one-stage discrete wavelet transform
(4) u2v.m: creates high-pass filter from low-pass filter
(5) v2u.m: creates low-pass filter from high-pass filter

This distribution also contains the data needed to recreate the results in 
the cited paper, and three MATLAB scripts that recreate the main results 
and demonstrate how to calculate data-driven wavelets with the above 
functions:
(1) exampleGaussianWhiteNoise.m: 
(2) exampleKS.m: 
(3) exampleHIT.m: 

As is, exampleHIT.m will take a long time to run. 

If you make use of this distribution, please cite "Discovering multiscale 
and self-similar structure with data-driven wavelets," by D. Floryan and 
M. D. Graham, PNAS, 2020. 
