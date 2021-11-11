function y = periodicNoise(len, f, blend, fs)
% Generates a pitched noise vector by looping a noise sequence.
% Parameters:
% len       : length of y (samples)
% f         : fundamental frequency of pitched noise
% blend     : periodic/white noise mix [0 (100% white) - 1 (100% periodic)]
% fs        : sample rate
%
% Author: Aaron Geldert
% Last modified: 22 October 2021

    % generate 1 period of noise
    T = round(fs/f);
    noise = rand(T,1) - 0.5;
    
    % repeat periodically
    numPeriods = ceil(len/T); 
    y = repmat(noise, numPeriods, 1);
    
    % 2nd order Low Pass Filter, 180 Hz cutoff
    Q = 0.7071;
    wc = 180 *2*pi/fs;
    alpha = 0.5 * sin(wc) / Q;
    beta = cos(wc);
    A = [1+alpha, -2*beta, 1-alpha];
    A = A/A(1);
    B = [(1-beta)/2, 1-beta, (1-beta)/2];
    B = B/A(1);
    y = filter(B,A,y);
    
    % truncate and apply gain
    y = blend * y(1:len,:);
    
    % blend in white noise and normalize
    y = y + (1-blend)*(rand(len,1)-0.5);
    y = y / max(abs(y));
end