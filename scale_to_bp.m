function [ScaledSignal] = scale_to_bp(signal, sbp, dbp)
%SCALE_TO_BP 
   % Function scales a signal of relative units into blood pressure units
   % according to the given maximum systolic pressure (sbp) and min
   % diastolic presure (dsp)
signal = signal-min(signal);
signal = signal/max(signal); 
ScaledSignal = signal *(sbp-dbp)+dbp;

   


end

