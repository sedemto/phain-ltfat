%% demo_FDGTusage: Demonstrating how to use "FDGT.m"
%   
%   This demo corresponds to "demo1_DGTusage.m" and illustrates the usage of
%   "FDGT.m" which is a faster but unreadable implementation of "DGT.m" [1].
%   See "demo1_DGTusage.m" and "DGT.m" for more details.
%   
%   The main difference between "DGT.m" and "FDGT.m" is generation of the
%   indices which is precomputed here so that repeated computation of the
%   same index is avoided. When computing DGT and invDGT to the same signal
%   repeatedly in a loop (such as in an optimization algorithm), FDGT and
%   invFDGT can save much computational time by calculating the indices
%   outside the loop (as demonstrated in this code).
%   
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)


%% Check MATLAB version and load signal

if verLessThan('matlab','9.1')
    error 'MATLAB 2016b or later is required due to repeated use of "implicit expansion"'
end
eval('load mtlb, signal = mtlb; clear mtlb') % loading speech signal


%% Set parameters of DGT (see "help DGT")

windowLen = 2^8;         % window length
shiftLen  = windowLen/4; % shifting stepsize
fftLen    = windowLen;   % number of FFT points

rotateFlag    = true; % deciding phase convention of DGT
zeroPhaseFlag = true; % deciding phase convention of window


%% Preprocess for DGT

currentFolder = cd; cd ../ % change folder for using these three functions

win = generalizedCosWin(windowLen,'blackman');   % generate window
dualWin = calcCanonicalDualWindow(win,shiftLen); % corresponding synthesis window

signal = zeroPaddingForDGT(signal,shiftLen,fftLen); % must run this before DGT

cd(currentFolder), clear currentFolder % come back to the current folder


%% Precompute indices to avoid unnecessary computation inside loop

[sigIdx,sumIdx,sumArray,ifftArray,rotIdx] = precomputationForFDGT(length(signal),windowLen,shiftLen,fftLen);

sigIdx = uint32(sigIdx); % this may slightly improve performance
sumIdx = uint32(sumIdx); % this may slightly improve performance
rotIdx = uint32(rotIdx); % this may slightly improve performance

if ~rotateFlag
    rotIdx = []; % existence of "rotIdx" indicates "rotateFlag" in FDGT
end


%% Calculate FDGT and its inverse repeatedly inside loop

result = signal;
maxIter = 100;

disp(['Running FDGT & invFDGT ' num2str(maxIter) ' times ...'])

tic
for n = 1:maxIter
    spec = FDGT(result,win,sigIdx,fftLen,rotIdx,zeroPhaseFlag);
    % do something ...
    result = invFDGT(spec,dualWin,sumIdx,sumArray,ifftArray,rotIdx,zeroPhaseFlag);
end
toc


%% Compare execution time with DGT and invDGT

currentFolder = cd; cd ../ % change folder for using these three functions

disp(['Running DGT & invDGT ' num2str(maxIter) ' times ...'])

tic
for n = 1:maxIter
    spec = DGT(result,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
    % do something ...
    result = invDGT(spec,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
end
toc

cd(currentFolder), clear currentFolder % come back to the current folder

