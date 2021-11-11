function out = astronaut(in,fs,model,type,k,sn,bn)
% Astronaut distortion audio effect
%   Usage: 
% out = astronaut(in,fs,model,type,k,sn,bn) 
% Parameters:
% in        : input data vector, fs
% fs    	: sample rate of in
% model 	: 1: Mercury, 2: Apollo
% type  	: nonlinearity type [1,2,3,4]
% k     	: nonlinearity coefficient, [0...1]
% sn    	: speech modulated noise gain, [0...1]
% bn    	: background noise gain, [0...1]
% out       : output data vector
%
% Example uses:
% mercury = astronaut(dry, fs, 1, 1, 0.1, 0.8, 0.2);
% apollo = astronaut(dry2, fs0, 2, 2, 0.6, 0.5, 0.15);
%
% Dependencies:
% - peakingFilter.m
% - periodicNoise.m
% 
% Authors: Tom Railio and Aaron Geldert
% Last Modified: 22 October 2021
% References: 
% [1]   S. Oksanen and V. Välimäki,
%   "Digital Modeling of the Vintage Telephone Sound", 
%   in Proc. ICMC2011, Huddersfield, UK, Jul. 31 - Aug. 5, 2011
% [2]   R. Bristow-Johnson, “Cookbook formulae for audio eq biquad filter 
%   coefficients,” http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt, 
%   2016, accessed: 2021-02-16.

%% Parameter Setup and Filter Design
if model <= 1
    % MERCURY MODELED SOUND
    % prepare peaking prefilters
    [B1,A1] = peakingFilter(290, -18, 0.7, fs);
    [B2,A2] = peakingFilter(1400, 18, 1.2, fs);
    % prepare peaking postfilter band
    [B3,A3] = peakingFilter(800, -12, 1.2, fs);
    postFilterFc = [600, 3800];
else
    % APOLLO MODELED SOUND
    % prepare peaking prefilters
    [B1,A1] = peakingFilter(90, -16, 0.8, fs);
    [B2,A2] = peakingFilter(950, 3, 1, fs);
    % prepare peaking postfilter band
    [B3,A3] = peakingFilter(3800, 11, 16, fs);
    postFilterFc = [55, 2200];
end

% POSTFILTER DESIGN: formulae from [2]
Q = 0.7071;
% 2nd order High Pass Filter
wc = postFilterFc(1) *2*pi/fs;
alpha = 0.5 * sin(wc) / Q;
beta = cos(wc);
A4 = [1+alpha, -2*beta, 1-alpha];
B4 = [(1+beta)/2, -1-beta, (1+beta)/2]/A4(1);
A4 = A4/A4(1);
% 2nd order Low Pass Filter
wc = postFilterFc(2) *2*pi/fs;
alpha = 0.5 * sin(wc) / Q;
beta = cos(wc);
A5 = [1+alpha, -2*beta, 1-alpha];
B5 = [(1-beta)/2, 1-beta, (1-beta)/2]/A5(1);
A5 = A5/A5(1);

%% Prefilter
% normalize signal
out = in./max(abs(in));

% apply 2 bands of peaking EQs
out = filter(B1,A1,out);
out = filter(B2,A2,out);

%% Nonlinearity
% Four potential waveshaping algorithms are implemented here
% Intensity parameter k : [0, 1]

switch type
    case 1
        % 2-sided sigmoid waveshaper (tube-like)
        k = k*50 + 1;
        out = 0.25/atan(k) * atan(k*out);
    case 2
        % 1-sided sigmoid waveshaper
        k = k*100 + 1;
        posInds = find(out>0.0);
        out(posInds) = 0.25/atan(k) * atan(k*out(posInds));
    case 3
        % 1-sided k-root waveshaper
        k = k*5 + 1;
        posInds = find(out>0.0);
        out(posInds) = out(posInds).^(1/k);
    case 4
        % carbon microphone waveshaper, see [1]
        if k > 0
            out = k^4*out.^5 + k^4*out.^4 - k^3*out.^4 - k^3*out.^3 +...
                k^2*out.^3 + k^2*out.^2 - k*out.^2 - k*out + out;
        end
end
% normalize
out = out./max(abs(out));

%% Speech Modulated Noise addition

% Envelope detector: uses a zero-phase LPF on rectified signal
[b,a] = butter(1,1400/fs); % 1st order LPF
env = filtfilt(b,a,abs(out));

% BP-filtered noise
[bb,aa] = butter(4,[500/(fs/2) 2700/(fs/2)]); % 4th order BPF
bp_noise = env.*filter(bb,aa,2*rand(length(in),1)-1); % apply filter and envelope
bp_noise = bp_noise./max(abs(bp_noise));  % normalize

% Noise applied
out = out + sn*(bp_noise);

%% Background Noise addition
if model <= 1
    % make Mercury modeled background noise
    humNoise = periodicNoise(length(out),814,0.95,fs);
else
    % make Apollo 11 modeled background hum
    humNoise = periodicNoise(length(out),57.5,0.96,fs);
end
% Noise applied
out = out + bn*humNoise;

%% Post-filtering
out = filter(B3,A3,out); % mid freq peaking band
out = filter(B4,A4,out); % HPF
out = filter(B5,A5,out); % LPF
out = filter(B5,A5,out); % LPF, applied twice for a 24 dB/octave slope

% normalize
out = out./max(abs(out));

end