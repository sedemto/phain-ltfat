function tightWin = calcCanonicalTightWindow(window,shiftLen)
% calcCanonicalTightWindow: Calculating the canonical tight window.
%   
%   Usage:
%      tightWin = calcCanonicalTightWindow(window,shiftLen);
%   
%   Input parameters:
%      window    : Input window (column vector).
%      shiftLen  : Shifting stepsize of DGT.
%   
%   Output parameters:
%      tightWin  : Associated canonical tight window.
%   
%   "tightWin = calcCanonicalTightWindow(window,shiftLen)" returns the canonical
%   tight window of the inputted window. This code is associated with "DGT.m"
%   implemented as a supporting material of [1]: It assumes that the number
%   of FFT points is not less than the length of window, "fftLen >= winLen".
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)

tightWin = buffer(window,shiftLen); % same as reshape (size: "shiftLen x something")
tightWin = tightWin./sqrt(sum(abs(tightWin).^2,2)); % inverting the frame operator
tightWin = tightWin(1:length(window)); % convert back into a vector
tightWin = tightWin(:); % ensure the output to be a column vector
end