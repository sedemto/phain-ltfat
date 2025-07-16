function spec = DGT(signal,window,shiftLen,fftLen,rotateFlag,zeroPhaseFlag)
% DGT: Index-based simple implementation of the discrete Gabor transform.
%   
%   Usage:
%      spec = DGT(signal,window,shiftLen);
%      spec = DGT(signal,window,shiftLen,fftLen);
%      spec = DGT(signal,window,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
%   
%   Input parameters:
%      signal        : Input signal (column vector).
%      window        : Analysis window (column vector).
%      shiftLen      : Shifting stepsize of DGT.
%      fftLen        : Number of FFT points.
%      rotateFlag    : Defining DGT (true: Eq. (7), false: Eq. (22) [1]).
%      zeroPhaseFlag : Defining window (true: zero-phased, false: linear-phased).
%   
%   Output parameters:
%      spec          : DGT coefficient (complex spectrogram).
%   
%   "spec = DGT(signal,window,shiftLen)" returns the discrete Gabor transform
%   of the inputted real-valued signal, where "window" is a real-valued
%   analysis window, and "shiftLen" is the shifting stepsize. By default,
%   the number of FFT points is the same as the length of window. It can be
%   adjusted by specifying "fftLen". As this code is merely an example for
%   demonstrating the time-frequency representation explained in [1], it is
%   assumed to be "fftLen >= winLen" for simplicity, where "winLen" is the
%   length of the window (note that it is possible to implement DGT in the 
%   case "fftLen < winLen" with a price of complicated inverse DGT).
%   
%   There are two flags for the index rotation explained in [1].
%   
%      "rotateFlag" activates the index rotation so that the phase is
%      represented as Eq. (7) in [1] (instead of Eq. (22)), i.e., the input
%      signal shares its time origin with the sinusoid.
%   
%      "zeroPhaseFlag" makes the window zero-phased, i.e., the linear-phase
%      component of the window is removed by the index rotation.
%   
%      By default, these flags are set to "false", which corresponds to DGT
%      of Eq. (22) in [1] (the window shares its time origin with the
%      sinusoid) with a window containing the linear-phase component.
%   
%   The signal must be modified by "zeroPaddingForDGT.m" before inputting to
%   this function for correct treatment of the periodic boundary condition [1].
%   
%   A faster implementation suitable for repeated use of DGT to the same
%   signal (as in optimization-based techniques) is available as "FDGT.m".
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)


%% Initial checks and default settings
narginchk(3,6)

sigLen = length(signal);
winLen = length(window);

if exist('fftLen','var') && ~isempty(fftLen)
    if fftLen < winLen, error 'Must satisfy "fftLen >= winLen"'; end
else
    fftLen = winLen;
end
if ~exist('rotateFlag','var') || isempty(rotateFlag)
    rotateFlag = false;
end
if ~exist('zeroPhaseFlag','var') || isempty(zeroPhaseFlag)
    zeroPhaseFlag = false;
end


%% Generating indices
% "idx" is the 2D index for cutting the signal into small time segments.
% "rotator" is the 2D index for rotation (by integer multiple of "shiftLen")
% within each segment.
[idx,rotator] = indexGenerator(sigLen,winLen,shiftLen,fftLen);


%% Windowing and zero-padding
% Zero-padding must be performed manually when the flags contain "true".
spec = [window.*signal(idx); zeros(fftLen-winLen,size(idx,2))];


%% Index rotation
% Periodic (or circular) shifting of the signal within each segment.
spec = reindexing(spec,winLen,zeroPhaseFlag,rotateFlag,rotator);


%% Calculating FFT and truncating negative frequency components
spec = fft(spec);
spec = spec(1:floor(fftLen/2)+1,:); % as the input signal is real-valued
end


%% Function for rotating signal within each segment
function spec = reindexing(spec,winLen,zeroPhaseFlag,rotateFlag,rotIdx)
if rotateFlag
    spec = spec(rotIdx); % convert representation from Eq. (22) to Eq. (7) of [1]
end
if zeroPhaseFlag
    spec = circshift(spec,-floor(winLen/2),1); % remove window's linear-phase component
end
end


%% Function for generating indices
function [sigIdx,rotIdx] = indexGenerator(sigLen,winLen,shiftLen,fftLen)
winIdx   = 0:winLen-1;          % row vector
winIdx   = winIdx(:);           % column vector
fftIdx   = 0:fftLen-1;          % row vector
fftIdx   = fftIdx(:);           % column vector
idxShift = 0:shiftLen:sigLen-1; % row vector

sigIdx = mod(winIdx + idxShift - floor(winLen/2), sigLen); % 2D index matrix
sigIdx = sigIdx + 1; % MATLAB-type index (starting from 1)

rotIdx = mod(fftIdx - idxShift, fftLen) + fftLen*(0:length(idxShift)-1);
rotIdx = rotIdx + 1; % MATLAB-type index (starting from 1)
end