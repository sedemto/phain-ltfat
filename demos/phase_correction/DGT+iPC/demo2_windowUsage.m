%% demo2_windowUsage: Demonstrating how to use window-related functions
%   
%   This demo illustrates how the window pairs work in terms of DGT and its
%   inverse transform [1]. It is recommended to read "demo1_DGTusage.m"
%   before running this code. Both canonical dual and tight windows are
%   utilized to see whether the specific window pair can reconstruct the
%   signal by DGT and invDGT. The reconstruction errors are displayed by
%   the texts in the command window.
%   
%   Figure 1: Canonical dual and tight window of the original window.
%      This figure shows the shapes of windows used in this code.
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
shiftLen  = windowLen/2; % shifting stepsize
fftLen    = windowLen*2; % number of FFT points

rotateFlag    = true; % deciding phase convention of DGT
zeroPhaseFlag = true; % deciding phase convention of window

signal = zeroPaddingForDGT(signal,shiftLen,fftLen); % must run this before DGT


%% Generate window and its canonical dual/tight version

win = generalizedCosWin(windowLen,'hann'); % original window

dualWin = calcCanonicalDualWindow(win,shiftLen);   % canonical dual window
tightWin = calcCanonicalTightWindow(win,shiftLen); % canonical tight window

figure, plot([win, dualWin, tightWin])
legend('original window', 'canonical dual window','canonical tight window')
codeOceanFigSave(gcf,'../../results','demo2_figure1.pdf') % save in Code Ocean


%% Original window cannot reconstruct signal (see "demo1_DGTusage.m")

spec = DGT(signal,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
reconst = invDGT(spec,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);

maxError = max(abs(signal-reconst))/max(abs(signal));
disp(['reconstruction error: ' num2str(maxError) ' (DGT: win, invDGT: win)'])


%% Dual window can reconstruct signal (see "demo1_DGTusage.m")

spec = DGT(signal,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
reconst = invDGT(spec,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);

maxError = max(abs(signal-reconst))/max(abs(signal));
disp(['reconstruction error: ' num2str(maxError) ' (DGT: win, invDGT: dualWin)'])


%% Dual window can reconstruct signal (see "demo1_DGTusage.m")

spec = DGT(signal,dualWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
reconst = invDGT(spec,win,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);

maxError = max(abs(signal-reconst))/max(abs(signal));
disp(['reconstruction error: ' num2str(maxError) ' (DGT: dualWin, invDGT: win)'])


%% Tight window can reconstruct signal (see "demo1_DGTusage.m")

spec = DGT(signal,tightWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);
reconst = invDGT(spec,tightWin,shiftLen,fftLen,rotateFlag,zeroPhaseFlag);

maxError = max(abs(signal-reconst))/max(abs(signal));
disp(['reconstruction error: ' num2str(maxError) ' (DGT: tightWin, invDGT: tightWin)'])

