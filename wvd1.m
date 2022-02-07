function [amb, tfrep] = wvd1(x,N);
%
% Wigner Ville Distribution
%
%  No windowing or time-resolution variablility
%
%
%
%
%
%
% Nathan Stevenson
% SPR June 2004


analytic_sig_ker = signal_kernal(x);
tfrep = real(1./N.*fft(ifftshift(analytic_sig_ker,1), N, 1));
amb = fftshift(1./N.*fft(analytic_sig_ker, N, 2),2);