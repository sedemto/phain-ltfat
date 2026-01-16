function IF = calcInstFreq(spec,diffSpec,fftLen,winLen,rotateFlag,flooringCoeff)
% calcInstFreq: Calculating instantaneous frequency at each bin.
%   
%   Usage:
%      IF = calcInstFreq(spec,diffSpec,fftLen,winLen);
%      IF = calcInstFreq(spec,diffSpec,fftLen,winLen,rotateFlag);
%      IF = calcInstFreq(spec,diffSpec,fftLen,winLen,rotateFlag,flooringCoeff);
%   
%   Input parameters:
%      spec          : Input DGT coefficient (complex spectrogram).
%      diffSpec      : DGT coefficient calculated by differentiated window.
%      fftLen        : Number of FFT points.
%      winLen        : Length of window used.
%      rotateFlag    : Flag used in "DGT.m" to compute "spec" and "diffSpec".
%      flooringCoeff : Small number for avoiding division by zero. 
%   
%   Output parameters:
%      IF            : Calculated bin-wise instantaneous frequency.
%   
%   "IF = calcInstFreq(spec,diffSpec,fftLen,winLen)" returns bin-wise relative
%   instantaneous frequency [1] calculated from the inputted spectrogram pair.
%   "diffSpec" must be obtained by the differentiated window corresponding to
%   the analysis window used for "spec". Such window pair can be picked up from
%   "generalizedCosWin.m" which includes several cosine-series-based windows
%   with their derivatives (e.g., Hann, Hamming, Blackman, and Nuttall). When
%   "rotateFlag" is false, the output becomes the absolute instantaneous
%   frequency instead of the relative instantaneous frequency [1].
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)

if ~exist('rotateFlag','var') || isempty(rotateFlag)
    rotateFlag = true;
end
if ~exist('flooringCoeff','var') || isempty(flooringCoeff)
    flooringCoeff = 1e-10;
end

powSpec = abs(spec).^2; % power spectrogram
powSpec = powSpec + flooringCoeff*max(powSpec(:)); % avoiding division by zero
% powSpec = max(powSpec,flooringCoeff*max(powSpec(:))); % this could be a choice

IF = -imag(diffSpec.*conj(spec)./powSpec); % calculating IF by Eq. (21) of [1]
IF = (fftLen/winLen)*IF; % compensation necessary when "fftLen ~= winLen"

if ~rotateFlag % wrapping effect must be treated separately if "rotateFlag == false"
    IF = (0:floor(fftLen/2))' + IF; % manually calculating Eq. (24) of [1]
end
end