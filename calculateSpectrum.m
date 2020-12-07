function [P1,f] = calculateSpectrum(Fourier,f_Hz)
%Calculates one sided frquency spectrum of fft transfomed Input Signal
    %Inputs are Fft transformed Signal as well as sampling frequency
    %Output is the one sided Spectrum and the corresponding frequency
    %vector

L=length(Fourier);
f=(0:(L/2))*f_Hz/L;
P1 = abs(Fourier(1:L/2+1) /L);
P1(2:end-1) = 2*P1(2:end-1);

end

