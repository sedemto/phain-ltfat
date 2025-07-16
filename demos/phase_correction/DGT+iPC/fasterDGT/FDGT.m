function spec = FDGT(signal,window,sigIdx,fftLen,rotIdx,zeroPhaseFlag)
% FDGT: Faster implementation of DGT.
%   
%   Usage:
%      spec = FDGT(signal,window,sigIdx);
%      spec = FDGT(signal,window,sigIdx,fftLen);
%      spec = FDGT(signal,window,sigIdx,fftLen,rotIdx,zeroPhaseFlag);
%   
%   Input parameters:
%      signal        : Input signal (column vector).
%      window        : Analysis window (column vector).
%      sigIdx        : Output of "precomputationForFDGT.m".
%      fftLen        : Number of FFT points.
%      rotIdx        : Output of "precomputationForFDGT.m".
%      zeroPhaseFlag : Defining window (true: zero-phased, false: linear-phased).
%   
%   Output parameters:
%      spec          : DGT coefficient (complex spectrogram).
%   
%   This is a faster but unreadable implementation of DGT [1] for repeated 
%   use to the same signal within a loop. See "DGT.m" for more details.
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)


%% Initial checks and default settings
narginchk(3,6)

winLen = length(window);

if exist('fftLen','var') && ~isempty(fftLen)
    if fftLen < winLen, error 'Must satisfy "fftLen >= winLen"'; end
else
    fftLen = winLen;
end
if ~exist('zeroPhaseFlag','var') || isempty(zeroPhaseFlag)
    zeroPhaseFlag = false;
end


%% Calculating DGT
if ~exist('rotIdx','var') || isempty(rotIdx)
    spec = window.*signal(sigIdx);
    if zeroPhaseFlag
        spec = [spec; zeros(fftLen-winLen,size(sigIdx,2))];
        spec = circshift(spec,-floor(winLen/2),1);
        spec = fft(spec);
    else
        spec = fft(spec,fftLen);
    end
else
    
    spec = [window.*signal(sigIdx); zeros(fftLen-winLen,size(sigIdx,2))];
    spec = spec(rotIdx);
    if zeroPhaseFlag
        spec = circshift(spec,-floor(winLen/2),1);
    end
    spec = fft(spec);
end

spec = spec(1:floor(fftLen/2)+1,:);
end