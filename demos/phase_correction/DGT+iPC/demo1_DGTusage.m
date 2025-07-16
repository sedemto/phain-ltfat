%% demo1_DGTusage: Demonstrating how to use "DGT.m"
%   
%   This demo illustrates the usage of "DGT.m" which is a supporting material
%   of [1]. The signal to be analyzed used in this demo is a female speech
%   installed in MATLAB by default (can be obtained by "load mtlb").
%   
%   Reconstruction error of DGT ("signal - invDGT(DGT(signal))") is displayed
%   by both figure and text in the command window.
%   
%   Figure 1: Spectrogram of the signal.
%      This figure shows the DGT coefficient of the signal (female speech).
%   
%   Figure 2: Original and reconstructed signal with their difference.
%      This figure shows time-domain signals with element-wise reconstruction
%      error. The original and reconstructed signal should coincide so that
%      only single signal is visible in the figure.
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
shiftLen  = windowLen/8; % shifting stepsize
fftLen    = windowLen;   % number of FFT points

rotateFlag    = true; % deciding phase convention of DGT
zeroPhaseFlag = true; % deciding phase convention of window


%% Generate window (see "help generalizedCosWin")

win = generalizedCosWin(windowLen,'blackman');


%% Preprocess for DGT

signal = zeroPaddingForDGT(signal,shiftLen,fftLen); % must run this before DGT


%% Calculate DGT

spec = DGT(signal,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);

figure, imagesc(20*log10(abs(spec))), colorbar, axis xy % display spectrogram
codeOceanFigSave(gcf,'../../results','demo1_figure1.pdf') % save in Code Ocean


%% Calculate inverse DGT

dualWin = calcCanonicalDualWindow(win,shiftLen); % corresponding synthesis window

reconst = invDGT(spec,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);


%% Evaluate reconstruction error "signal - invDGT(DGT(signal))"

figure, plot([signal, reconst, signal-reconst])
legend('original signal', 'reconstructed signal','reconstruction error')
codeOceanFigSave(gcf,'../../results','demo1_figure2.pdf') % save in Code Ocean

maxError = max(abs(signal-reconst))/max(abs(signal));
disp(['relative maximum reconstruction error: ' num2str(maxError)])

