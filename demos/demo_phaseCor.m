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
%% LTFAT DGT setup
% setup tight window
g = gabtight(wtype, a, M, w);
if wtype == "hann"
    % derivative of Hann window
    x = (0:w-1)'/(w);
    g_diff = -0.5*sin(2*pi.*x)*max(g);
else 
    % for other window functions
    g_diff = numericalDiffWin(g);
end

% setup DGT and invDGT  using LTFAT
G_ltfat = @(x) comp_sepdgtreal(x, g, a, M, phasetype);
G_adj_ltfat = @(u) comp_isepdgtreal(u, g, size(u,2)*a, a, M, phasetype);
G_diff = @(x) comp_sepdgtreal(x, g_diff, a, M, phasetype);

% Compute the DGT of the segment with a gap
sgram = G_ltfat(segment.gapped);

%% calculate instantaneous frequency
omega = @(x) calcInstFreq(G_ltfat(x), G_diff(x), M, w);
omega_y = omega(segment.gapped);
%% phase correction original
R = @(z, omega) instPhaseCorrection(z, omega, a, M);
R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);

%% phase correction new
phaseRotations = @(omega) rotations(omega, a, M);
phaseCor = exp(-1i*phaseRotations(omega_y));
invPhaseCor = exp(1i*phaseRotations(omega_y));
R_new = @(z) phaseCor.*z;
R_adj_new =  @(z) invPhaseCor.*z;

%% Test
numberOfTests = 10;
numberOfIterations =10000;
test_times = table('Size',[numberOfTests 9],'VariableTypes',{'double','double','double','double','double','double','double','double','double'});
test_times.Properties.VariableNames = ["original code (only R)","LTFAT code (only R)","R improvement [%]","original code (only R*)","LTFAT code (only R*)","R* improvement [%]","original code (both R and R*)","LTFAT code (both R and R*)","improvement [%]"];
for num=1:numberOfTests
    % testing R
    tic
    for i=1:numberOfIterations
        a = R(sgram,omega_y);
    end
    test_times{num,1} = toc;
    tic
    for i=1:numberOfIterations
        b = R_new(sgram);
    end
    test_times{num,2} = toc;
    
    % testing R*
    tic
    for i=1:numberOfIterations
        c = R_adj(sgram,omega_y);
    end
    test_times{num,4} = toc;
    tic
    for i=1:numberOfIterations
        d = R_adj_new(sgram);
    end
    test_times{num,5} = toc;
    
    % testing both R and R*
    tic
    for i=1:numberOfIterations
        e_cor = R(sgram,omega_y);
        e_invCor = R_adj(sgram,omega_y);
    end
    test_times{num,7} = toc;
    tic
    for i=1:numberOfIterations
        f_cor = R_new(sgram);
        f_invCor = R_adj_new(sgram);
    end
    test_times{num,8} = toc;
end
test_times{:,3} = 100 - test_times{:,2}./test_times{:,1}*100;
test_times{:,6} = 100 - test_times{:,5}./test_times{:,4}*100;
test_times{:,9} = 100 -test_times{:,8}./test_times{:,7}*100;

disp(test_times)

%% functions
% MIT License
% 
% Copyright (c) 2019 Kohei Yatabe
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)
%
%   The below functions were adapted from [1]

function IF = calcInstFreq(spec,diffSpec,fftLen,winLen,flooringCoeff)
    % calc_IF: Calculating instantaneous frequency at each bin.
    
    if ~exist('flooringCoeff','var') || isempty(flooringCoeff)
        flooringCoeff = 1e-10;
    end
    
    powSpec = abs(spec).^2; % power spectrogram
    powSpec = powSpec + flooringCoeff*max(powSpec(:)); % avoiding division by zero
    
    IF = -imag(diffSpec.*conj(spec)./powSpec); % calculating IF by Eq. (21) of [1]
    IF = (fftLen/winLen)*IF; % compensation necessary when "fftLen ~= winLen"
end

function iPCspec = instPhaseCorrection(spec,IF,shiftLen,fftLen)
    % instPhaseCorrection: Calculating instantaneous-phase-corrected spectrogram.
    
    sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
    freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize
    
    idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
    cumPhase = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
    
    iPCspec = exp(-1i*cumPhase).*spec; % implementation of Eq. (29) of [1]
end

function spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen)
    % invInstPhaseCorrection: Inverting instantaneous phase correction.
    
    sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
    freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize
    
    idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
    cumPhase = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
    
    spec = exp(1i*cumPhase).*iPCspec; % inverting phase correction
end

function phaseCor = rotations(IF,shiftLen,fftLen)
    % Calculating phase rotations for phase-corrected spectrogram.
    
    sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
    freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize
    
    idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
    phaseCor = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
=======
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
%% LTFAT DGT setup
% setup tight window
g = gabtight(wtype, a, M, w);
if wtype == "hann"
    % derivative of Hann window
    x = (0:w-1)'/(w);
    g_diff = -0.5*sin(2*pi.*x)*max(g);
else 
    % for other window functions
    g_diff = numericalDiffWin(g);
end

% setup DGT and invDGT  using LTFAT
G_ltfat = @(x) comp_sepdgtreal(x, g, a, M, phasetype);
G_adj_ltfat = @(u) comp_isepdgtreal(u, g, size(u,2)*a, a, M, phasetype);
G_diff = @(x) comp_sepdgtreal(x, g_diff, a, M, phasetype);

% Compute the DGT of the segment with a gap
sgram = G_ltfat(segment.gapped);

%% calculate instantaneous frequency
omega = @(x) calcInstFreq(G_ltfat(x), G_diff(x), M, w);
omega_y = omega(segment.gapped);
%% phase correction original
R = @(z, omega) instPhaseCorrection(z, omega, a, M);
R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);

%% phase correction new
phaseRotations = @(omega) rotations(omega, a, M);
phaseCor = exp(-1i*phaseRotations(omega_y));
invPhaseCor = exp(1i*phaseRotations(omega_y));
R_new = @(z) phaseCor.*z;
R_adj_new =  @(z) invPhaseCor.*z;

%% Test
numberOfTests = 10;
numberOfIterations =10000;
test_times = table('Size',[numberOfTests 9],'VariableTypes',{'double','double','double','double','double','double','double','double','double'});
test_times.Properties.VariableNames = ["original code (only R)","LTFAT code (only R)","R improvement [%]","original code (only R*)","LTFAT code (only R*)","R* improvement [%]","original code (both R and R*)","LTFAT code (both R and R*)","improvement [%]"];
for num=1:numberOfTests
    % testing R
    tic
    for i=1:numberOfIterations
        a = R(sgram,omega_y);
    end
    test_times{num,1} = toc;
    tic
    for i=1:numberOfIterations
        b = R_new(sgram);
    end
    test_times{num,2} = toc;
    
    % testing R*
    tic
    for i=1:numberOfIterations
        c = R_adj(sgram,omega_y);
    end
    test_times{num,4} = toc;
    tic
    for i=1:numberOfIterations
        d = R_adj_new(sgram);
    end
    test_times{num,5} = toc;
    
    % testing both R and R*
    tic
    for i=1:numberOfIterations
        e_cor = R(sgram,omega_y);
        e_invCor = R_adj(sgram,omega_y);
    end
    test_times{num,7} = toc;
    tic
    for i=1:numberOfIterations
        f_cor = R_new(sgram);
        f_invCor = R_adj_new(sgram);
    end
    test_times{num,8} = toc;
end
test_times{:,3} = 100 - test_times{:,2}./test_times{:,1}*100;
test_times{:,6} = 100 - test_times{:,5}./test_times{:,4}*100;
test_times{:,9} = 100 -test_times{:,8}./test_times{:,7}*100;

disp(test_times)

%% functions
% MIT License
% 
% Copyright (c) 2019 Kohei Yatabe
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)
%
%   The below functions were adapted from [1]

function IF = calcInstFreq(spec,diffSpec,fftLen,winLen,flooringCoeff)
    % calc_IF: Calculating instantaneous frequency at each bin.
    
    if ~exist('flooringCoeff','var') || isempty(flooringCoeff)
        flooringCoeff = 1e-10;
    end
    
    powSpec = abs(spec).^2; % power spectrogram
    powSpec = powSpec + flooringCoeff*max(powSpec(:)); % avoiding division by zero
    
    IF = -imag(diffSpec.*conj(spec)./powSpec); % calculating IF by Eq. (21) of [1]
    IF = (fftLen/winLen)*IF; % compensation necessary when "fftLen ~= winLen"
end

function iPCspec = instPhaseCorrection(spec,IF,shiftLen,fftLen)
    % instPhaseCorrection: Calculating instantaneous-phase-corrected spectrogram.
    
    sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
    freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize
    
    idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
    cumPhase = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
    
    iPCspec = exp(-1i*cumPhase).*spec; % implementation of Eq. (29) of [1]
end

function spec = invInstPhaseCorrection(iPCspec,IF,shiftLen,fftLen)
    % invInstPhaseCorrection: Inverting instantaneous phase correction.
    
    sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
    freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize
    
    idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
    cumPhase = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
    
    spec = exp(1i*cumPhase).*iPCspec; % inverting phase correction
end

function phaseCor = rotations(IF,shiftLen,fftLen)
    % Calculating phase rotations for phase-corrected spectrogram.
    
    sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
    freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize
    
    idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
    phaseCor = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
>>>>>>> 13786de3beae7075cc56dfb6ac19ec5323a91963
end