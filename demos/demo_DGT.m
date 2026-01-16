close all
clear
clc
ltfatstart
rng(0)
addpath("../U-PHAIN/")
addpath(genpath("phase_correction"))
%% load signal
path = "../dataset/DPAI_originals/audio_original_example0.wav";
[data,fs] = audioread(path);

%% parameters

% parameter settings for STFT/DGT
w = 2048; % window length
a = w/4; % hop size
M = w; % number of freq. rows
wtype = 'hann'; % window type
phasetype = 0; % 0-freqinv, 1-timeinv

%% create gap and select segment around it
% usually only a small segment around a gap is selected
% to simulate it a gap needs to be made
% this is a simplyfied version of the original code

gapLength = 50; % [ms]
h = round(fs*gapLength/1000); % [samples]

notDegraded = 0.5; % [s]
segment.length = round((length(data) - 2*notDegraded*fs));
idx = round(notDegraded*fs) + (1:segment.length);
% making a gap at random     
s = round((w + 1) + rand()*(segment.length - 2*w - h));
f = s + h - 1;
segment.mask = true(segment.length, 1);
segment.mask(s:f) = false;
segment.data = data(idx);
segment.gapped = segment.data.*segment.mask;

% shortening the segment
[firstIdx, L] = shortenForDGT(w, a, s, f);
origL = L;
enoughL = ceil(L/lcm(a, M))*lcm(a, M);
if L < enoughL
    L = enoughL;
end
 
lastIdx = firstIdx + L - 1;
if firstIdx >= 1 && lastIdx <= segment.length
    segment.mask = segment.mask(firstIdx:lastIdx);
    segment.gapped = segment.gapped(firstIdx:lastIdx);
    segment.data = segment.data(firstIdx:lastIdx);
    idx = idx(firstIdx:firstIdx + L - 1);
    st = 1;
    ed = L;
else
    firstIdx = max(firstIdx, 1);
    padding = zeros(lastIdx - length(segment.data), 1);
    segment.mask = [segment.mask(firstIdx:end); true(size(padding))];
    segment.gapped = [segment.gapped(firstIdx:end); padding];
    segment.data = [segment.data(firstIdx:end); padding];
    idx = idx(firstIdx:firstIdx + length(segment.data) - length(padding) - 1);
    st = 1;
    ed = L - length(padding);
end
%% LTFAT DGT setup
% setup tight window
g = gabtight(wtype, a, M, w);

% setup DGT and invDGT  using LTFAT
G_ltfat = @(x) comp_sepdgtreal(x, g, a, M, phasetype);
G_adj_ltfat = @(u) comp_isepdgtreal(u, g, size(u,2)*a, a, M, phasetype);

%% original DGT setup

zeroPhaseFlag = true;
rotateFlag = true;

% setup tight window
[win, ~] = generalizedCosWin(w, 'hanning');
tight_win = calcCanonicalTightWindow(win, a);
tight_win = tight_win/norm(tight_win)*sqrt(a/w);

% setup DGT and invDGT original code
[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(size(segment.gapped,1), w, a, M);

G = @(x) FDGT(x, tight_win, sigIdx, M, rotIdx, zeroPhaseFlag);
G_adj = @(u) invFDGT(u, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*w;



%% Test
numberOfTests = 10;
numberOfIterations =10000;
test_times = table('Size',[numberOfTests 9],'VariableTypes',{'double','double','double','double','double','double','double','double','double'});
test_times.Properties.VariableNames = ["original code (only DGT)","LTFAT code (only DGT)","DGT improvement [%]","original code (only inverse DGT)","LTFAT code (only inverse DGT)","invDGT improvement [%]","original code (both DGT and invDGT)","LTFAT code (both DGT and invDGT)","improvement [%]"];
for num=1:numberOfTests
    % testing DGT
    tic
    for i=1:numberOfIterations
        a = G(segment.gapped);
    end
    test_times{num,1} = toc;
    tic
    for i=1:numberOfIterations
        b = G_ltfat(segment.gapped);
    end
    test_times{num,2} = toc;
    
    % testing inverse DGT
    temp1 = G(segment.gapped);
    temp2 = G_ltfat(segment.gapped);
    tic
    for i=1:numberOfIterations
        c = G_adj(temp1);
    end
    test_times{num,4} = toc;
    tic
    for i=1:numberOfIterations
        d = G_adj_ltfat(temp2);
    end
    test_times{num,5} = toc;
    
    % testing both DGT and inverse DGT
    tic
    for i=1:numberOfIterations
        e = G(segment.gapped);
        e_sig = G_adj(e);
    end
    test_times{num,7} = toc;
    tic
    for i=1:numberOfIterations
        f = G_ltfat(segment.gapped);
        f_sig = G_adj_ltfat(f);
    end
    test_times{num,8} = toc;
end
test_times{:,3} = 100 - test_times{:,2}./test_times{:,1}*100;
test_times{:,6} = 100 - test_times{:,5}./test_times{:,4}*100;
test_times{:,9} = 100 -test_times{:,8}./test_times{:,7}*100;

disp(test_times)
