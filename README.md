# DDWD
Data-driven wavelets

README for Data-Driven Wavelet Decomposition (DDWD)

This distribution contains the code needed to create data-driven wavelets, 
as described in "Discovering multiscale and self-similar structure with 
data-driven wavelets," by D. Floryan and M. D. Graham, PNAS, 2020. 

This distribution contains five primary MATLAB functions: 
* waveletOpt.m: calculates data-driven wavelets
* dwtos.m: one-stage discrete wavelet transform
* idwtos.m: inverse one-stage discrete wavelet transform
* u2v.m: creates high-pass filter from low-pass filter
* v2u.m: creates low-pass filter from high-pass filter

This distribution also contains the data needed to recreate the results in 
the cited paper, and three MATLAB scripts that recreate the main results 
and demonstrate how to calculate data-driven wavelets with the above 
functions:
* exampleGaussianWhiteNoise.m: recreates main results for Gaussian white noise data
* exampleKS.m: recreates main results for Kuramoto-Sivashinsky data
* exampleHIT.m: recreates main results for turbulence data

As is, exampleHIT.m will take a long time to run. 

If you make use of this distribution, please cite "Discovering multiscale 
and self-similar structure with data-driven wavelets," by D. Floryan and 
M. D. Graham, PNAS, 2020. 
