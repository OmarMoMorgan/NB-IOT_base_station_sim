function Hd = test
%TEST Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.11 and DSP System Toolbox 9.13.
% Generated on: 30-Dec-2023 08:07:18

% FIR Window Lowpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 7680000;  % Sampling Frequency

N    = 100;      % Order
Fc   = 200000;   % Cutoff Frequency
flag = 'scale';  % Sampling Flag

% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
