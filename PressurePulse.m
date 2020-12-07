%% Lena Castel-Wohnlich and Wolfgang Fuchs present: 
% Determination of ejection duration and flow model using pressure curves

%% Load Data, define important variables
clear;
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

x1=filtfilt(IIR_Lowpass(f_c_L).sosMatrix, IIR_Lowpass(f_c_L).ScaleValues,struct(1).pressure);
x2=filtfilt(IIR_Lowpass(f_c_L).sosMatrix, IIR_Lowpass(f_c_L).ScaleValues,struct(2).pressure);
x3=filtfilt(IIR_Lowpass(f_c_L).sosMatrix, IIR_Lowpass(f_c_L).ScaleValues,struct(3).pressure);
figure; hold on
subplot(3,1,1);plot(t,struct(1).pressure,t,x1)
subplot(3,1,2);plot(t,struct(2).pressure,t,x2)
subplot(3,1,3);plot(t,struct(3).pressure,t,x3)


