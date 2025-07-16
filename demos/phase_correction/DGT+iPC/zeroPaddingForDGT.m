function sigLong = zeroPaddingForDGT(signal,shiftLen,fftLen)
% zeroPaddingForDGT: Zero-padding required for treating DGT with periodic boundary.
%   
%   Usage:
%      sigLong = zeroPaddingForDGT(signal,shiftLen,fftLen);
%   
%   Input parameters:
%      signal   : Input signal (column vector).
%      shiftLen : Shifting stepsize of DGT.
%      fftLen   : Number of FFT points.
%   
%   Output parameters:
%      sigLong  : Zero-padded signal.
%   
%   "sigLong = zeroPaddingForDGT(signal,shiftLen,fftLen)" returns the signal
%   with zero-padding required for treating DGT with periodic boundary [1].
%   This function must be run before using "DGT.m" for avoiding the error in
%   the index generator (local function inside "DGT.m").
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)

sigLen = length(signal);
c = lcm(shiftLen,fftLen); % least common multiple
zeroPadLen = ceil(sigLen/c)*c - sigLen; % length of zeros for padding
sigLong = [signal; zeros(zeroPadLen,1)]; % zero-padding
end