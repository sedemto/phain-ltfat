%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   PHASE-AWARE AUDIO INPAINTING (PHAIN)                  %
%           remade using the LTFAT library (only includes U-PHAIN)        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT License
% 
% Copyright (c) 2023 TomoroTanaka
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

%   [1] Tanaka, Tomoro, Kohei Yatabe, and Yasuhiro Oikawa,"PHAIN: Audio
%       inpainting via phase-aware optimization with instantaneous frequency,"
%       IEEE/ACM Transactions on Audio, Speech, and Language Processing, Sep 2024.
%
% all sections except "functions" were adapted from main.m in [1]
%
% the main difference between LTFAT PHAIN and the original PHAIN code is in
% the "parameters" section, specifically the setup of a tight window and DGT

close all
clear
clc
ltfatstart
rng(0)

addpath(genpath('../dataset/'))
%% loading
soundDir = "../dataset/DPAI_originals/";
ext = ".wav";
Sounds = dir(soundDir + "*" + ext);
NN = length(Sounds);
data = cell(NN,1);
info = audioinfo(Sounds(1).name);
fs = info.SampleRate;
for nn = 1:NN
    data{nn} = audioread(Sounds(nn).name);
end
clear audio info

%% settings
     
gaps = 5:5:50; % [ms]
gapNum = length(gaps);

N = 8; % # of gaps

methodLabels = {'U_PHAIN'};

for i = 1:length(methodLabels)
    solution.(methodLabels{i}) = cell(NN, gapNum);  % initialization of restored results
end

SNR  = NaN(NN, gapNum, N, length(methodLabels));  % SNR per gap
TIME = NaN(NN, gapNum, N, length(methodLabels));  % execution time per gap
SNR_procedure = cell(NN, gapNum, N, length(methodLabels));  % iterative behavior

%% parameters

% parameter settings for STFT/DGT
w = 2048; % window length
a = w/4; % hop size
M = w; % number of freq. rows
wtype = 'hann'; % window type
phasetype = 0; % 0-freqinv, 1-timeinv

% update input of the Chambolle-Pock alg.
param.updateInputCP = false; % default false (in the original code) 

% setup tight window and its derivative 
g = gabtight(wtype, a, M, w);

if wtype == "hann"
    % derivative of Hann window
    x = (0:w-1)'/(w);
    g_diff = -0.5*sin(2*pi.*x)*max(g);
else 
    % for other window functions
    g_diff = numericalDiffWin(g);
end
% g_diff = numericalDiffWin(g);

% setup DGT, invDGT and derivative of DGT using LTFAT
param.G = @(x) comp_sepdgtreal(x, g, a, M, phasetype);
param.G_adj = @(u) comp_isepdgtreal(u, g, size(u,2)*a, a, M, phasetype);
param.G_diff = @(x) comp_sepdgtreal(x, g_diff, a, M, phasetype);

% definition of instantaneous frequency (omega)
param.omega = @(x) calcInstFreq(param.G(x), param.G_diff(x), M, w);

% def.of phase correction (R) and time-directional difference (D)
param.R = @(z, omega) instPhaseCorrection(z, omega, a, M);
param.R_adj =  @(z, omega) invInstPhaseCorrection(z, omega, a, M);
param.D = @(z) z(:,1:end-1) - z(:,2:end);
param.D_adj = @(z) [z(:,1), (z(:,2:end) - z(:,1:end-1)), -z(:,end)];


% settings for generalized CP algorithm
paramsolver.epsilon = 0.01;  % for stopping criterion

paramsolver.tau = 0.25;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.alpha = 1;  % relaxation parameter

paramsolver.lambda = 1; % threshold (regularization parameter)
%% inpainting

for nn = 1:NN

    signal = data{nn};

    for m = 1:gapNum

        fprintf('\nSignal: %d / %d', nn, NN)
        fprintf('\nGap Length: %d [ms]\n', gaps(m))

        gapLength = gaps(m); % [ms]
        h = round(fs*gapLength/1000); % [samples]
        full.length = length(signal);
        full.mask = true(full.length, 1);

        notDegraded = 0.5; % [s]
        segment.length = round((length(signal) - 2*notDegraded*fs) / N);

        for i = 1:length(methodLabels)
            solution.(methodLabels{i}){nn, m} = signal;
        end

        for n = 1:N
            
            fprintf('\nGap Number: %d / %d\n', n, N)
            idx = round(notDegraded*fs) + ((n - 1)*segment.length+1:n*segment.length);

            % making a gap      
            s = round((w + 1) + rand()*(segment.length - 2*w - h));
            f = s + h - 1;
            segment.mask = true(segment.length, 1);
            segment.mask(s:f) = false;
            full.mask(idx) = segment.mask;
           
            segment.data = signal(idx);
            segment.max = max(abs(segment.data));
            segment.data = segment.data/segment.max;
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%% U-PHAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('U-PHAIN...\n')
            tic
            param.type = 'U';
            paramsolver.I = 100;
            paramsolver.J = 10;

            [segment.solution, SNR_procedure{nn, m, n, 1}] = ...
                UPHAIN_ltfat(segment.gapped, segment.mask, param, paramsolver, segment.data);
            
                    
            solution.U_PHAIN{nn, m}(idx) = segment.solution(st:ed)*segment.max;
            TIME(nn, m, n, 1) = toc;

        end
   
        % calculating SNR
        fprintf('\nevaluation start\n')
        for i = 1:length(methodLabels)
            restored = solution.(methodLabels{i}){nn, m};
            groundtruth = signal(~full.mask);
            result = restored(~full.mask);
            for n = 1:N
                SNR(nn, m, n, i) =  snr(groundtruth(1 + h*(n - 1):h*n),groundtruth(1 + h*(n - 1):h*n)-result(1 + h*(n - 1):h*n));
            end
        end
        fprintf('\nevaluation done!\n')

    end

end

%% plot

snr_vec = squeeze(median(SNR, [1,3]));
figure(Position = [614 157 873 830])
p = plot(snr_vec, LineWidth = 2);
grid on
legend({"U-PHAIN"}, Interpreter = 'latex')
xlabel('gap length 5:5:50 [ms]', Interpreter = 'latex')
xticks = 1:10;
xticklabels = {"5", "10", "15", "20", "25", "30", "35", "40", "45", "50"};
ylabel('SNR at gaps [dB]', Interpreter = 'latex')
ax = gca;
ax.FontSize = 15;
ax.TickLabelInterpreter = 'latex';


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

function diffWin = numericalDiffWin(window,zeroPadLen)
    % Zero-padding for alleviating the periodic boundary effect
    if ~exist('zeroPadLen','var') || isempty(zeroPadLen)
        zeroPadLen = 0;
    end
    longWin = [window; zeros(zeroPadLen,1)]; % zero-padding
    
    % Generating the frequency vector (index) for spectral derivative
    winLen = length(window);
    longLen = length(longWin);
    M = floor((longLen-1)/2);
    
    fftIdx = ifftshift([zeros(mod(longLen-1,2)),-M:M]); % index generation
    fftIdx = fftIdx(:)*winLen/longLen; % normalization
    
    % Calculating spectral derivative
    diffWin = ifft(1i*fftIdx.*fft(longWin),'symmetric'); % spectral method
    diffWin = diffWin(1:winLen); % truncating padded zeros
end
