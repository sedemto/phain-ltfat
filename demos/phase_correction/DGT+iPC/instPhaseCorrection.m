function iPCspec = instPhaseCorrection(spec,IF,shiftLen,fftLen)
% instPhaseCorrection: Calculating instantaneous-phase-corrected spectrogram.
%   
%   Usage:
%      iPCspec = instPhaseCorrection(spec,IF,shiftLen,fftLen);
%   
%   Input parameters:
%      spec     : Input DGT coefficient (complex spectrogram).
%      IF       : Bin-wise estimate of instantaneous frequency.
%      shiftLen : Shifting stepsize of DGT.
%      fftLen   : Number of FFT points.
%   
%   Output parameters:
%      iPCspec  : Instantaneous-phase-corrected (iPC) spectrogram.
%   
%   "iPCspec = instPhaseCorrection(spec,IF,shiftLen,fftLen)" returns the
%   instantaneous-phase-corrected (iPC) spectrogram [1]. The input complex
%   spectrogram "spec" must be calculated by "DGT.m" with "rotateFlag" being
%   "true" (DGT must be defined as Eq. (7) of [1] because this code implements
%   iPC in Eq. (29) of [1]). The instantaneous frequency (IF) can be estimated
%   by using "calcInstFreq.m" which directly obtains IF based on its standard
%   definition. Other estimator of IF can be utilized to improve performance.
%   Some applications of iPC spectrogram are listed in [1].
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)

sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize

idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
cumPhase = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value

iPCspec = exp(-1i*cumPhase).*spec; % implementation of Eq. (29) of [1]
end