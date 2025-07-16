close all
clear
clc
ltfatstart
rng(0)
addpath("../U-PHAIN/")
%% load signal
path = "../dataset/DPAI_originals/audio_original_example0.wav";
[data,fs] = audioread(path);

%% parameters

% parameter settings for STFT/DGT
w = 2048; % window length
a = w/4; % hop size
M = w; % number of freq. rows
wtype = 'hann'; % window type
gapLength = 50; % [ms]

%% create gap and select segment around it
% usually only a small segment around a gap is selected
% to simulate it a gap needs to be made
% this is a simplyfied version of the original code

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
%% Test
corrupted = segment.gapped;
mask = segment.mask;

proj_original_func = @(x) projGamma(x, mask, corrupted);
proj_new_func = @(x) projGamma_new(x, mask, corrupted);
proj_new = @(x) x.*(1-mask) + corrupted.*mask;

% projection using original projection function
tic
for i=1:10000
    p1 = proj_original_func(corrupted);
end
toc

% projection using new projection function (fair comparison)
tic
for i=1:10000
    p2 = proj_new_func(corrupted);
end
toc

% projection without function
tic
for i=1:10000
    p3 = proj_new(corrupted);
end
toc


%% Functions

function x = projGamma(x,mask , corrupted)
    x(mask) = corrupted(mask);
end

function x = projGamma_new(x,mask , corrupted)
    x = x.*(1-mask) + corrupted.*mask;
end