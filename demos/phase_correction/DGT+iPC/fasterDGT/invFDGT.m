function signal = invFDGT(spec,window,sumIdx,sumArray,ifftArray,rotIdx,zeroPhaseFlag)
% invFDGT: Faster implementation of invDGT.
%   
%   Usage:
%      signal = invFDGT(spec,window,sumIdx,sumArray,ifftArray);
%      signal = invFDGT(spec,window,sumIdx,sumArray,ifftArray,rotIdx,zeroPhaseFlag);
%   
%   Input parameters:
%      spec          : Input DGT coefficient (complex spectrogram).
%      window        : Synthesis window (column vector).
%      sumIdx        : Output of "precomputationForFDGT.m".
%      sumArray      : Output of "precomputationForFDGT.m".
%      ifftArray     : Output of "precomputationForFDGT.m".
%      rotIdx        : Output of "precomputationForFDGT.m".
%      zeroPhaseFlag : Defining window (true: zero-phased, false: linear-phased).
%   
%   Output parameters:
%      signal        : Reconstructed signal.
%   
%   This is a faster but unreadable implementation of invDGT [1] for repeated 
%   use to the same signal within a loop. See "invDGT.m" for more details.
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)

%   Author: Kohei Yatabe (2019)


%% Initial checks and default settings
narginchk(5,7)

winLen = length(window); % assuming the lengths of dual-window pair are the same

if ~exist('zeroPhaseFlag','var') || isempty(zeroPhaseFlag)
    zeroPhaseFlag = false;
end


%% Calculating inverse DGT
ifftArray(1:size(spec,1),:) = spec;
sig = ifft(ifftArray,'symmetric');

if zeroPhaseFlag
    sig = circshift(sig,floor(winLen/2),1);
end
if exist('rotIdx','var') && ~isempty(rotIdx) 
    sig(rotIdx) = sig;
end

sig = window.*sig(1:winLen,:); % multiplying synthesis window

sumArray(sumIdx) = sig(:); % resolve overlapping
signal = sum(sumArray,2);  % overlap add
end