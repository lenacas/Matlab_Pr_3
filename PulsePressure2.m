%% Lena Castel-Wohnlich and Wolfgang Fuchs present: 
% Determination of ejection duration and flow model using pressure curves

%% 0.Load Data, define important variables
clear;
close all;
load("Project3_02_Data.mat");
fs=250;                             %Sampling f; all freq in Hz
L= length(struct(1).pressure);      %Luckily all measurments of equal L
t= (0:L-1) *1/fs;

%% 1.Filter high f noise
% Spectrum suggests high frequency noise above 50 Hz.
% Initial jR filter approach dropped due to high signal distortion and
% subsequent problems with accuratly splitting the signal and detecting. 
% Might be intresting for later comparison though?

% FIR Lowpass Filter designed with filtertoolbox
temp=FIR_Lowpass;
for i=1:3
struct(i).filtered_signal= filtfilt(temp.Numerator,1,struct(i).pressure);
end

%% 2.Split signal into beats
% using findpeaks, we detect the lowest values, where we assume a single
% beat starts.

% Split Filtered Signal
for i = 1:3
    [peak, location_f] = findpeaks(-struct(i).filtered_signal,t,'MinPeakDistance', 0.45);
    location_f = location_f *fs/1;
   
for j = 1:(length(location_f)-1)
    %filtered single beats creation and saved to struct
    struct(i).f_s_beat(j).signal = struct(i).filtered_signal(round(location_f(j)):round(location_f(j+1))); % round to create integer
    struct(i).f_s_beat(j).time = (0:length(struct(i).f_s_beat(j).signal)-1) *1/fs;
end

% Split Unfiltered Signal (task 8.)
    [peak, location] = findpeaks(-struct(i).pressure,t,'MinPeakDistance', 0.5);
    location = location *fs/1;

    if length(location_f) == length(location)
        for j = 1:(length(location)-1)
    %filtered single beats creation by cutting at lowest values  and saved to struct
    struct(i).s_beat(j).signal = struct(i).pressure(round(location(j)):round(location(j+1))); % round to create integer
    struct(i).s_beat(j).time = (0:length(struct(i).s_beat(j).signal)-1) *1/fs;
        end
    end
    if length(location_f) ~= length(location)
        for j = 1:(length(location)-2)
    %filtered single beats creation by cutting at lowest values  and saved to struct
    struct(i).s_beat(j).signal = struct(i).pressure(round(location(j)):round(location(j+1))); % round to create integer
    struct(i).s_beat(j).time = (0:length(struct(i).s_beat(j).signal)-1) *1/fs;
        end 
    end
end 


%% 3. Scale signal
% by function: subtract minimum, divide by max, multiply with delta bp, 
% add diastolic bp;
for i=1:3
    for j = 1:(length(struct(i).f_s_beat))  
      struct(i).f_s_beat(j).signal= scale_to_bp(struct(i).f_s_beat(j).signal,struct(i).sbp,struct(i).dbp);     
    end
end

%% 4. Three-point-moving-average filter 
% ASK! Order of assignment keept, but wouldnt it be smarter to filter and
% subsequently scale in order to avoid value distortion by the filter
F=[1,1,1]/3;
for i=1:3
    for j = 1:(length(struct(i).f_s_beat)-1)
       struct(i).f_s_beat(j).signal = filter(F,1,struct(i).f_s_beat(j).signal- struct(i).f_s_beat(j).signal(1))+struct(i).f_s_beat(j).signal(1);
    end    
end

% was passiert mit 1. und letzten Datenpunkt? Zeitliche Verschiebung


%% 5.+ 6. Find and save Systolic peak and Dicrotic notch
% Systolic peak = abs Maximum
% Dicrotic notch equivalent to first local minimum after systolic peak

% time fängt bei 0 an, loc bei 1 (t= 0.004 s)

peakprominence=3; 
for i = 1:3
    temp=0;
for j = 1:length(struct(i).f_s_beat)
    [temp_max,f_s_locmax(i,j)] = max(struct(i).f_s_beat(j).signal);
    [temp_dic,f_s_ejtimes(i,j)] = findpeaks(-struct(i).f_s_beat(j).signal,struct(i).f_s_beat(j).time,'MinPeakProminence',peakprominence);
    temp= temp+f_s_ejtimes(i,j);
end
end

f_s_locmax = f_s_locmax*1/fs;

%% 7. Ejection Time
% in seconds
for i = 1:3
    struct(i).ejtimeav_f = 0;
for j = 1:length(struct(i).f_s_beat)
    struct(i).ejtimeav_f = struct(i).ejtimeav_f + f_s_ejtimes(i,j);
end
    struct(i).ejtimeav_f = struct(i).ejtimeav_f/j;
end

%% 8. Repeat steps for unfiltered Signal
% Splitting of unfiltered signal is done in 2) (more structured)
for i=1:3
    struct(i).ejtimeav=0;
    for j = 1:(length(struct(i).s_beat))
        struct(i).s_beat(j).signal= scale_to_bp(struct(i).s_beat(j).signal,struct(i).sbp,struct(i).dbp);
        struct(i).s_beat(j).signal= filter(F,1,struct(i).s_beat(j).signal- struct(i).s_beat(j).signal(1))+struct(i).s_beat(j).signal(1);
       [temp_max,s_locmax(i,j)] = max(struct(i).s_beat(j).signal);
       [temp,s_ejtimes(i,j)] = findpeaks(-struct(i).s_beat(j).signal,struct(i).s_beat(j).time,'MinPeakProminence',peakprominence);
       struct(i).ejtimeav = struct(i).ejtimeav + s_ejtimes(i,j);
    end
    struct(i).ejtimeav = struct(i).ejtimeav/j;
end

s_locmax = s_locmax*1/fs;

%% 9.Bland Altman Plots
% Compare each beat for each patient - average required
% 

for i=1:3
        figure;
        BlandAltman(s_locmax(i,1:length(struct(i).f_s_beat)),f_s_locmax(i,1:length(struct(i).f_s_beat)),2,i,'Location of Systolic Peak');
        figure;
        BlandAltman(s_ejtimes(i,1:length(struct(i).f_s_beat)),f_s_ejtimes(i,1:length(struct(i).f_s_beat)),2,i,'Ejection Time');
end 

% The differences are small. The mean of differences are for all three
% patients in a range of milli seconds or even smaller (systolic peak of #2
% and ejection time of #3).
% All data differences are in the confidence intervall of 95%.
% It shows, that the used filter keeps the reliability of the original
% data. 
% For a better comparison, more data points should be used here.

%% 10.Average filtered single beats

%determine max index of the single beat vectors per ID
max_length =[0 0 0];
for i=1:3
for j=1:(length(struct(i).f_s_beat))
    if length(struct(i).f_s_beat(j).signal)>max_length(i)
        max_length(i) =length(struct(i).f_s_beat(j).signal);
    end    
end
end

for i=1:3 
    % initialize Nan matrix of correct Dimension to fit all the beats
    temp = NaN(length(struct(i).f_s_beat),max_length(i)) ;
    
for j=1:(length(struct(i).f_s_beat))
    % fill Nan Matrix with values
    temp(j,1:length(struct(i).f_s_beat(j).signal))=struct(i).f_s_beat(j).signal;
end
    % averages signals and saves them to struct (NaN ignored)
    struct(i).mean=mean(temp,'omitnan');
end

%% 11.Inflection point (first zero crossing of 2nd Derivative)

for i=1:3
    %find values were 2nd Derivative is negative
    temp=find(diff(diff(struct(i).mean))<0);
    InflectionPoint(i)=temp(1);
    %first negative value assumed as deflection point index, else subtract 
    %one and take the first positive value
    
    % verschiebung?? Welche Datenpunkte gehen verloren?? Wie muss man
    % verschieben?
end

%% 12. Plot Mean beat signal and Triangle
for i=1:3
    figure; hold on; 
    yyaxis left;
    plot((0:length(struct(i).mean)-1)/fs ,struct(i).mean)
    xlabel("Time [s]")
    ylabel("Blood Pressure [mmHg]")
    title("Triangular flow model ID #"+num2str(i))
    subtitle("of "+num2str(length(struct(i).f_s_beat))+" averaged hearbeats")
    yyaxis right;
    ylabel("Flow [%]")
    line([0 (InflectionPoint(i)-1)/fs],[0 100])
    line ([(InflectionPoint(i)-1)/fs struct(i).ejtimeav_f],[100 0])
    legend("Averaged Signal","Pressure Flow")
end