%% TASK 2: Prodosic Features

% PATHS NEEDED
addpath ./packages/voicebox/
addpath ./packages/SpeechProsodyAnalysisTools/SpeechProsodyAnalysisTools/F0_HarmonicAnalysis/


%% LOAD DATA and downsample
load lab2_data.mat;

% resample
Fs = 48000; % Original sampling frequency 48kHz
Fs_down = 11025; % sampling frequency for the mfcc extraction
y = resample(speech_sample,Fs_down,Fs); % downsample

%% Calculate Short Time Energy
% Calculate energy contour
ste = conv(y.^2,hamming(0.01*Fs_down)); % short time energy using a hamming window of 0.01s

% clip the contour (half a frame from the beginning and half a frame from the end)
ste = ste(0.01*Fs_down/2: end-0.01*Fs_down/2)

% Plot the speech sample waveform and the STE contour
t_ste = linspace(0,length(ste),5); % time vector for plotting the original signal and the STE 
t = linspace(0,2.5*Fs,6); % time vector for plotting the waveform
subplot(2,1,1),plot(speech_sample)
xlabel('Time(s)')
ylabel('Amplitude')
axis([0 2.5*Fs -1 1])
set(gca,'xtick',t)
set(gca,'xticklabel',get(gca,'xtick')/Fs)
title('Raw speech signal')

subplot(2,1,2),plot(ste)
xlabel('Time(s)')
ylabel('Energy')
set(gca,'xticklabel',get(gca,'xtick')/10^4)
title('Short Term Energy (STE) of  the speech signal')
% Extract STE features according to the instruction
mean_ste = mean(ste)
std_ste = std(ste)
prctile10_ste = prctile(ste,10)
prctile90_ste = prctile(ste,90)
kurtosis_ste = kurtosis(ste)
%% Extract pitch contour
[F0,strength,T_ind,wflag] = getF0(y, Fs_down);% extract F0 contour, include unvoiced/silent data

% Plot the raw speech sample waveform and the pitch contour
t_F0 = linspace(0,length(F0),5);% time vector for plotting the F0_contour
t = linspace(0,2.5*Fs,6); % time vector for plotting the waveform
subplot(2,1,1),plot(speech_sample)
xlabel('Time(s)')
ylabel('Amplitude')
set(gca,'xtick',t)
set(gca,'xticklabel',get(gca,'xtick')/Fs)
axis([0 2.5*Fs -1 1])
title('Raw speech signal')

subplot(2,1,2),plot(F0)
xlabel('Time(s)')
ylabel('Frequency/(Hz)')
set(gca,'xticklabel',get(gca,'xtick')/100)
title('F0-contour of  the speech signal')

% Extract 5 F0 features according to the instruction
voiced_F0 = F0(F0>0);
mean_F0 = mean(voiced_F0)
std_F0 = std(voiced_F0)
prctile10_F0 = prctile(voiced_F0,10)
prctile90_F0 = prctile(voiced_F0,90)
kurtosis_F0 = kurtosis(voiced_F0)

%% Extract durations 
% Try to design alogrithm by yourself, using the hint given in the
% instruction
F0 = F0>0;
seg = diff(F0);
frame = find(seg~=0);
voiceEnds = [];
unvoiceBegs = [];
voiceBegs = [];
unvoiceEnds = [];
for i = frame
    if seg(i) == -1
        voiceEnds = [voiceEnds, i];
        unvoiceBegs = [unvoiceBegs, i+1];
    elseif seg(i) == 1
        unvoiceEnds = [unvoiceEnds, i];
        voiceBegs = [voiceBegs, i+1];
    end
end
if F0(1) == 1
    voiceBegs = [1, voiceBegs];
elseif F0(1) == 0
    unvoiceBegs = [1, unvoiceBegs];
end

if F0(end) == 1
    voiceEnds = [voiceEnds, length(F0)];
elseif F0(end) == 0
    unvoiceEnds = [unvoiceEnds, length(F0)];
end
length_voiced = (voiceEnds - voiceBegs +1) *0.01;
length_unvoiced = (unvoiceEnds - unvoiceBegs +1)*0.01;

% Extract durations and voiced ratio according to the instruction
mean_voiced = mean(length_voiced)
std_voiced = std(length_voiced)

mean_unvoiced = mean(length_unvoiced)
std_unvoiced = std(length_unvoiced)
voicing_ratio = sum(length_voiced)*100/length(F0)

%% Check task 2 results
% check that the calcualted features correspond with the first row of the provided data matric (precalculated feature data)
% Your extracted feature, 5 from F0, 5 from STE, and 5 from Durations correspondingly
proso_features = [mean_F0 std_F0 prctile10_F0 prctile90_F0 kurtosis_F0 mean_ste std_ste prctile10_ste prctile90_ste kurtosis_ste mean_voiced std_voiced mean_unvoiced std_unvoiced voicing_ratio]
precalculated_proso = training_data_proso(1,:) 