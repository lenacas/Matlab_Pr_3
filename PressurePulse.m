%% Lena Castel-Wohnlich and Wolfgang Fuchs present: 
% Determination of ejection duration and flow model using pressure curves

%% Load Data, define important variables
clear;
close all;
load("Project3_02_Data.mat");
fs=250;                             %Sampling f; all freq in Hz
L= length(struct(1).pressure);      %Luckily all measurments of equal L


%% Visualize data
% Plot data
figure; hold on 
t= (0:L-1) *1/fs;
plot(t,struct(1).pressure,t,struct(2).pressure,t,struct(3).pressure)
legend("P1","P2","P3")

% Plot spectra
[P1,f1]=calculateSpectrum(fft(struct(1).pressure),fs);
[P2,f2]=calculateSpectrum(fft(struct(2).pressure),fs);
[P3,f3]=calculateSpectrum(fft(struct(3).pressure),fs);
figure; hold on
subplot(3,1,1); semilogy(f1,P1); title("FFT of Signal #1"); xlabel("Frequency [Hz]")
subplot(3,1,2); semilogy(f2,P2); title("FFT of Signal #2"); xlabel("Frequency [Hz]")
subplot(3,1,3); semilogy(f3,P3); title("FFT of Signal #3"); xlabel("Frequency [Hz]")

%% Filter 
% FFT suggests signal below and noise above 20Hz; 60Hz peak= powerline noise; 
% What is the 40Hz peak ?? Signal or noise? 40Hz noise?
% Filter toolbox: Creation of Lowpass with fc =30
% -> IIR filter with minimum order selcted (eliptic window)
f_c_L=30;

struct(1).x=filtfilt(IIR_Lowpass(f_c_L).sosMatrix, IIR_Lowpass(f_c_L).ScaleValues,struct(1).pressure);
struct(2).x=filtfilt(IIR_Lowpass(f_c_L).sosMatrix, IIR_Lowpass(f_c_L).ScaleValues,struct(2).pressure);
struct(3).x=filtfilt(IIR_Lowpass(f_c_L).sosMatrix, IIR_Lowpass(f_c_L).ScaleValues,struct(3).pressure);
figure; hold on
subplot(3,1,1);plot(t,struct(1).pressure,t,struct(1).x); xlabel("time (s)")
subplot(3,1,2);plot(t,struct(2).pressure,t,struct(2).x); xlabel("time (s)")
subplot(3,1,3);plot(t,struct(3).pressure,t,struct(3).x); xlabel("time (s)")

%% Split to single beats and scale plots
% using findpeaks, we detect the lowest values, where we assume a single
% beat starts. 
% 1 Pa = 133,322 mmHg

for i = 1:3
    
[struct(i).peak,struct(i).location] = findpeaks(-struct(i).x,t,'MinPeakDistance', 0.5);
struct(i).location = struct(i).location *fs/1;

for ii = 1:(length(struct(i).location)-1)
    struct(i).singlebeat(ii).signal = struct(i).x(struct(i).location(ii):struct(i).location(ii+1));
    struct(i).singlebeat(ii).time = (0:length(struct(i).singlebeat(ii).signal)-1) *1/fs;
end
end 

mmHg_Pa = 0.00750062; % change of unit --> I don't know which unit it is


figure; title("Single beats #1");hold on;
for ii = 1:length(struct(1).location)-1
    plot(struct(1).singlebeat(ii).time,struct(1).singlebeat(ii).signal*mmHg_Pa); xlabel("time (s)")
    % ylim([struct(1).dbp struct(1).sbp]) --> gggggrrrr, why does it not
    % work?
end


figure; title("Single beats #2");
for ii = 1:length(struct(2).location)-1
    hold on; 
    plot(struct(2).singlebeat(ii).time,struct(2).singlebeat(ii).signal*mmHg_Pa); xlabel("time (s)")
    % ylim([struct(2).dbp struct(2).sbp])
end

figure; title("Single beats #3");
for ii = 1:length(struct(3).location)-1
    hold on; 
    plot(struct(3).singlebeat(ii).time,struct(3).singlebeat(ii).signal*mmHg_Pa); xlabel("time (s)")
    % ylim([struct(3).dbp struct(3).sbp])
end



% figure; hold on
% subplot(3,1,1); plot(struct(1).x); hold on; plot(struct(1).location,-struct(1).peak,'x')
% subplot(3,1,2); plot(struct(2).x); hold on; plot(struct(2).location,-struct(2).peak,'x')
% subplot(3,1,3); plot(struct(3).x); hold on; plot(struct(3).location,-struct(3).peak,'x')
% xlabel("time (s)")
% this was just for testing if it fits

