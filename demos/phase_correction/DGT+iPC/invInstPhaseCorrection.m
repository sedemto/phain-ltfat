function spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen)
% invInstPhaseCorrection: Inverting instantaneous phase correction.
%   
%   Usage:
%      spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen);
%   
%   Input parameters:
%      iPCspec  : Instantaneous-phase-corrected (iPC) spectrogram.
%      IF       : Bin-wise estimate of instantaneous frequency.
%      shiftLen : Shifting stepsize of DGT.
%      fftLen   : Number of FFT points.
%   
%   Output parameters:
%      spec     : Complex spectrogram (without phase correction).
%   
%   "spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen)" returns complex
%   spectrogram corresponding to the inputted instantaneous-phase-corrected
%   (iPC) spectrogram by inverting the phase correction [1]. All parameters
%   must be the same as those used for "instPhaseCorrection.m" to recover the
%   original spectrogram.
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)

sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize

idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
cumPhase = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
%size(exp(1i*cumPhase))
%size(iPCspec)
spec = exp(1i*cumPhase).*iPCspec; % inverting phase correction
end