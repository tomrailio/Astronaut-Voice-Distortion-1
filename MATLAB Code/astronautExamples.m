%% Astronaut Voice Distortion Effect Examples
% This script loads an audio sample, processes it using the astronaut
% function, and plays the output audio. A spectrogram plot is also
% included as an example.
%
% Created by: Aaron Geldert and Tom Railio
% Last modified: 22 October 2021

% Load an audio signal
[x, fs] = audioread('step_dry.wav');

%% Example: MERCURY 6
y1 = astronaut(x, fs, 1, 1, 0.1, 0.8, 0.2);
sound(y1,fs);
% audiowrite("mercury_demo.wav",y1,fs); % optional save to .wav file

%%  Example: APOLLO 11
y2 = astronaut(x, fs, 2, 1, 0.3, 0.7, 0.2);
sound(y2,fs);
% audiowrite("apollo_demo.wav",y2,fs); % optional save to .wav file

%% Spectrogram plot
subplot(311);
spectrogram(x,hamming(2048),1024,8092,fs,'yaxis')
ylim([20/1000 4.1])
title('Input Signal Spectrogram');

subplot(312);
spectrogram(y1,hamming(2048),1024,8092,fs,'yaxis')
ylim([20/1000 4.1])
title('Mercury model of Astronaut Voice Distortion');

subplot(313);
spectrogram(y2,hamming(2048),1024,8092,fs,'yaxis')
ylim([20/1000 4.1])
colormap(jet)
title('Apollo model of Astronaut Voice Distortion');
