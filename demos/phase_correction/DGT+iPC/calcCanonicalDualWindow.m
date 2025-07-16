function dualWin = calcCanonicalDualWindow(window,shiftLen)
% calcCanonicalDualWindow: Calculating the canonical dual window.
%   
%   Usage:
%      dualWin = calcCanonicalDualWindow(window,shiftLen);
%   
%   Input parameters:
%      window   : Input window (column vector).
%      shiftLen : Shifting stepsize of DGT.
%   
%   Output parameters:
%      dualWin  : Associated canonical dual window.
%   
%   "dualWin = calcCanonicalDualWindow(window,shiftLen)" returns the canonical
%   dual window of the inputted window. This code is associated with "DGT.m"
%   implemented as a supporting material of [1]: It assumes that the number
%   of FFT points is not less than the length of window, "fftLen >= winLen".
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)

dualWin = buffer(window,shiftLen); % same as reshape (size: "shiftLen x something")
dualWin = dualWin./sum(abs(dualWin).^2,2); % inverting the frame operator
dualWin = dualWin(1:length(window)); % convert back into a vector
dualWin = dualWin(:); % ensure the output to be a column vector
end