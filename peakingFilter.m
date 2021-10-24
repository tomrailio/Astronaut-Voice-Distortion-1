function [B, A] = peakingFilter(fc, dBgain, Q, fs)
% Calculates Biquad coefficients for a 2nd order peaking filter.
% Parameters:
% fc        : center frequency, Hz
% dBgain    : peak gain, decibels
% Q         : quality factor
% fs        : sample rate
%
% Reference: RBJ Audio Eq Cookbook
% https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
%
% Author: Aaron Geldert
% Last modified: 9 February 2021

    gain = sqrt(10^(dBgain/20));
    wc = 2*pi*fc/fs;
    c = cos(wc);
    s = sin(wc);
    alpha = 0.5*s/Q;
    
    a0 = 1 + alpha/gain;
    a1 = -2*c;
    a2 = 1 - alpha/gain;
    b0 = 1 + alpha*gain;
    b1 = a1;
    b2 = 1 - alpha*gain;
    
    A = [a0 a1 a2] / a0;
    B = [b0 b1 b2] / a0;
end