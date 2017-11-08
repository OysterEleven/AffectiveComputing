%% TASK 1: MFCC

% PATHS NEEDED
addpath ./packages/voicebox/
addpath ./packages/SpeechProsodyAnalysisTools/SpeechProsodyAnalysisTools/F0_HarmonicAnalysis/


%% LOAD DATA and downsample
load lab2_data.mat;
whos

% resample
Fs = 48000; % Original sampling frequency 48kHz
Fs_down = 11025; % sampling frequency for the mfcc extraction
y = resample(speech_sample,Fs_down,Fs); % downsample


%% Apply Pre-emphasize Filter
% pre-emphasis filter
% Define the numerator and denominator coefficients for the rational transfer function.
% Please refer to the instruction and reference page of filter (command: help filter)
prefilter_b = [1,-0.98]; % filter coefficients b, which is the numerator
prefilter_a = [1]; % filter coefficients, which is the denominator
 
prefiltered_sample = filter(prefilter_b,prefilter_a,y); % filtering

% Plot the pre-emphasis filter frequency response in Hz using the filter visualization tool ?vtool?
%%%%%fvtool(prefilter_b, prefilter_a, 'Fs', Fs_down) % pre-emphasis filter visualization

 
%% Extract MFCC Feature
mfcc_c = melcepst(prefiltered_sample, Fs_down); % MFCC extraction


%% Plot Your MFCC Countours and Waveforms
len1 = length(speech_sample)
t = linspace(0,2.5*Fs,6); % time vector for plotting the waveform
t_mfcc = linspace(0,2.5); % time vector for plotting mfcc values
subplot(2,1,1),plot(speech_sample)
xlabel('Time(s)')
ylabel('Amplitude')
axis([0 2.5*Fs -1 1])
set(gca,'xtick',t)
set(gca,'xticklabel',get(gca,'xtick')/Fs)
title('Raw speech signal')

subplot(2,1,2),plot(mfcc_c)
xlabel('Time(s)')
ylabel('Amplitude')
set(gca,'xticklabel',get(gca,'xtick')/100)
title('MFCC contours of the speech signal')

%% Calculate the 'mean' of each MFCC feature contour to produce a total of 12 features
% check that the calcualted features correspond with the first row of the provided data matric 
%(precalculated feature data)
mean_mfcc = mean(mfcc_c)
precalculated_mfcc = training_data_mfcc(1,:)
difference = mean_mfcc - precalculated_mfcc