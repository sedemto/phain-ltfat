function signal = invDGT(spec,window,shiftLen,fftLen,rotateFlag,zeroPhaseFlag)
% invDGT: Index-based simple implementation of the inverse discrete Gabor transform.
%   
%   Usage:
%      signal = invDGT(spec,window,shiftLen);
%      signal = invDGT(spec,window,shiftLen,fftLen);
%      signal = invDGT(spec,window,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
%   
%   Input parameters:
%      spec          : Input DGT coefficient (complex spectrogram).
%      window        : Synthesis window (column vector).
%      shiftLen      : Shifting stepsize of DGT.
%      fftLen        : Number of FFT points.
%      rotateFlag    : Defining DGT (true: Eq. (7), false: Eq. (22) [1]).
%      zeroPhaseFlag : Defining window (true: zero-phased, false: linear-phased).
%   
%   Output parameters:
%      signal        : Reconstructed signal.
%   
%   "signal = invDGT(spec,window,shiftLen)" returns the inverse discrete Gabor
%   transform of the complex spectrogram obtained by "DGT.m", where "window"
%   is the corresponding synthesis window, and "shiftLen" is the shifting
%   stepsize. This code assumes that the length of the synthesis window is
%   the same as that of the analysis window for avoiding bothersome calculations.
%   All inputted parameters (except the window) are assumed to be the same
%   as those used in "DGT.m" to obtain "spec". The synthesis window must be
%   a dual window of the analysis window for perfect reconstruction [1].
%   
%   A faster implementation suitable for repeated use of invDGT to the same
%   signal (as in optimization-based techniques) is available as "invFDGT.m".
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)


%% Initial checks and default settings
narginchk(3,6)

winLen = length(window); % assuming the lengths of dual-window pair are the same
sigLen = size(spec,2)*shiftLen;

if ~exist('fftLen','var') || isempty(fftLen)
    if floor(winLen/2)+1 ~= size(spec,1), error 'Please specify "fftLen"', end
    fftLen = winLen; % "spec" is assumed to be obtained by "DGT.m" without "fftLen"
end
if ~exist('rotateFlag','var') || isempty(rotateFlag)
    rotateFlag = false;
end
if ~exist('zeroPhaseFlag','var') || isempty(zeroPhaseFlag)
    zeroPhaseFlag = false;
end


%% Generating indices
% "idx" is the 2D index relating the time segment into the original signal.
% "rotator" is the 2D index for canceling the rotation imposed by "DGT.m".
% "idx2" is the 2D index for handling overlapping parts of the signal.
[idx,rotator,idx2] = indexGenerator(sigLen,winLen,shiftLen,fftLen,size(spec,2));


%% Calculating inverse FFT
% "symmetric" treats the spectrum conjugate-symmetric so that the converted
% signal is ensured to be real-valued.
signal = ifft([spec; zeros(fftLen-size(spec,1),size(spec,2))],'symmetric');


%% Inverting index rotation
signal = reindexingBackward(signal,winLen,zeroPhaseFlag,rotateFlag,rotator);


%% Multiplying synthesis window
signal = window.*signal(1:winLen,:);


%% Overlap add (in shortcut form)
% "sparse" creates a sparse matrix whose column has the same dimension as
% the original signal. Its row handles overlapping parts which are added by "sum".
signal = full(sum(sparse(idx(:),idx2(:),signal(:)),2));
end


%% Function for inverting the rotation imposed by DGT
function signal = reindexingBackward(signal,winLen,zeroPhaseFlag,rotateFlag,rotIdx)
if zeroPhaseFlag
    signal = circshift(signal,floor(winLen/2),1); % add window's linear-phase component
end
if rotateFlag
    signal(rotIdx) = signal; % invert representation from Eq. (7) to Eq. (22) [1]
end
end


%% Function for generating indices
function [sigIdx,rotIdx,sumIdx] = indexGenerator(sigLen,winLen,shiftLen,fftLen,timeLen)
winIdx   = 0:winLen-1;          % row vector
winIdx   = winIdx(:);           % column vector
fftIdx   = 0:fftLen-1;          % row vector
fftIdx   = fftIdx(:);           % column vector
idxShift = 0:shiftLen:sigLen-1; % row vector

sigIdx = mod(winIdx + idxShift - floor(winLen/2), sigLen); % 2D index matrix
sigIdx = sigIdx + 1; % MATLAB-type index (starting from 1)

rotIdx = mod(fftIdx - idxShift, fftLen) + fftLen*(0:length(idxShift)-1);
rotIdx = rotIdx + 1; % MATLAB-type index (starting from 1)

sumIdx = 0:timeLen-1;
sumIdx = sumIdx + 1; % MATLAB-type index (starting from 1)
sumIdx = repmat(sumIdx,winLen,1);
end