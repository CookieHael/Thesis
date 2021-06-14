function [DCT_cost, WHT_cost, Haar_wavelet_cost] = digital_transformations_costs(N)
%DIGITAL_TRANSFORMATIONS_COSTS Outputs an array with the number of
%operations for some common transformations, based on size of transform
DCT_op = 17/9 * N * log2(N) - 17/27*N - 1/9 * (-1)^(log2(N))*log2(N) + 7/54 * (-1)^(log2(N)) + 3/2;  %flops  N^2 coeff (size?)  DCT2  https://doi.org/10.1016/j.sigpro.2008.01.004
WHT_add = N*log(N); %additions only, N^2 (1 bit) coeff
Haar_wavelet_mult = N; %2N-2 additions, N^2 coeff  Fundamentals of the Discrete Haar Wavelet Transform Duraisamy Sundararajan
Haar_wavelet_add = 2*N-2;

DCT_cost = 1* DCT_op + 1 * N^2 * 8;
WHT_cost = 1*WHT_add + 1 * N^2;
Haar_wavelet_cost = 1 * Haar_wavelet_mult + 1 * Haar_wavelet_add + 1 * N^2 * 2;

end

