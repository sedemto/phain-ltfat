function diffWin = numericalDiffWin(window,zeroPadLen)
% numericalDiffWin: Numerically differentiating a window.
%   
%   Usage:
%      diffWin = numericalDiffWin(window);
%      diffWin = numericalDiffWin(window,zeroPadLen);
%   
%   Input parameters:
%      window     : Input window (column vector).
%      zeroPadLen : Length of zeros for zero-padding.
%   
%   Output parameters:
%      diffWin    : Numerically differentiated window.
%   
%   "diffWin = numericalDiffWin(window)" returns a numerically differentiated
%   window calculated by the spectral method [1]. As it imposes the periodic
%   boundary condition, zero-padding (activated by specifying "zeroPadLen")
%   is implemented to alleviate the boundary effect.
%   
%   This function is included for calculating the instantaneous frequency
%   using "calcInstFreq.m" [2] with a discretely defined window. For a
%   generalized cosine window (including Hann, Hamming, Blackman, Nuttall, etc.),
%   use "generalizedCosWin.m" which implements analytic derivatives.
%   
%   [1] Lloyd N. Trefethen, Spectral Methods in MATLAB, SIAM, 2000.
%   
%   [2] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)


%% Zero-padding for alleviating the periodic boundary effect
if ~exist('zeroPadLen','var') || isempty(zeroPadLen)
    zeroPadLen = 0;
end
longWin = [window; zeros(zeroPadLen,1)]; % zero-padding


%% Generating the frequency vector (index) for spectral derivative
winLen = length(window);
longLen = length(longWin);
M = floor((longLen-1)/2);

fftIdx = ifftshift([zeros(mod(longLen-1,2)),-M:M]); % index generation
fftIdx = fftIdx(:)*winLen/longLen; % normalization


%% Calculating spectral derivative
diffWin = ifft(1i*fftIdx.*fft(longWin),'symmetric'); % spectral method
diffWin = diffWin(1:winLen); % truncating padded zeros
end