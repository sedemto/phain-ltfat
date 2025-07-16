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
% some sections adapted from main.m in [1]
%
% the main difference between LTFAT PHAIN and the original PHAIN code is in
% the "parameters" section, specifically the setup of a tight window and DGT

close all
clear
clc
ltfatstart
rng(0)

addpath("../U-PHAIN")
addpath(genpath('../dataset/'))
addpath(genpath("phase_correction"))
%% loading
soundDir = "../dataset/DPAI_originals/*0";
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
     
gaps = 30; % [ms]
gapNum = length(gaps);

N = 1; % # of gaps

methodLabels = {'U_PHAIN'};

for i = 1:length(methodLabels)
    solution.(methodLabels{i}) = cell(NN, gapNum);  % initialization of restored results
end

SNR  = NaN(NN, gapNum, N, length(methodLabels));  % SNR per gap
SNR_procedure = cell(NN, gapNum, N, length(methodLabels));  % iterative behavior
SNR_procedure_PHAIN = cell(NN, gapNum, N, length(methodLabels));  % iterative behavior
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

g_diff = numericalDiffWin(g);

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

            param.type = 'U';
            paramsolver.I = 100;
            paramsolver.J = 10;
            disp("Time of the LTFAT PHAIN")
            tic
            % for i=1:100
            [segment.solution, SNR_procedure{nn, m, n, 1}] = ...
                UPHAIN_ltfat(segment.gapped, segment.mask, param, paramsolver, segment.data);
            % end
            toc
            
            param.a = a;
            param.w = w;
            param.M = M;
            disp("Time of the original PHAIN")
            tic
            % for i=1:100
            [segment.solution_PHAIN, SNR_procedure_PHAIN{nn, m, n, 1}] = ...
                PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);
            % end
            toc
            % plot only the gap
            figure;
            plot(segment.gapped(~segment.mask), 'Color',[0 0 0 0]); hold on;
            plot(segment.solution(~segment.mask)); hold on;
            plot(segment.solution_PHAIN(~segment.mask)); hold on;
            plot(segment.data(~segment.mask)); hold on;
            legend("","LTFAT PHAIN","original PHAIN","groundtruth")
            xlabel("Samples in gap")
            ylabel("Amplitude")
            
            solution.U_PHAIN{nn, m}(idx) = segment.solution(st:ed)*segment.max;

            snr_original = snr(segment.data(~segment.mask),segment.data(~segment.mask)-segment.solution_PHAIN(~segment.mask));
            snr_ltfat = snr(segment.data(~segment.mask),segment.data(~segment.mask)-segment.solution(~segment.mask));

            disp("SNR of the reconstruction the original code: " + snr_original)
            disp("SNR of the reconstruction using the LTFAT code: " + snr_ltfat)
        end
    end
end
