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
function [outsig, snr_procedure] = PHAINmain(insig, mask, param, paramsolver, oracle)

    % param
    %   .a
    %   .M
    %   .w
    %   .type 
    
    
    % paramsolver
    %   .sigma .... step size
    %   .tau ...... step size
    %   .alpha .... relaxation paramter
    %   .lambda ... threshold
    %   .epsilon .. stop criterion
    %   .x0 ....... initial value of primal variable
    %   .u0 ....... initial value of dual variable
    %   .I ........ number of inner iterations
    %   .J ........ number of outer iterations
    
    
    %% iPC DGT
    
    a = param.a;
    M = param.M;
    w = param.w;
    
    [win, ~] = generalizedCosWin(w, 'hanning');
    tight_win = calcCanonicalTightWindow(win, a);
    tight_win = tight_win/norm(tight_win)*sqrt(a/w);
    diff_win = numericalDiffWin(tight_win);
        
    zeroPhaseFlag = true;
    rotateFlag = true;
    
    [sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(length(insig), w, a, M);
    
    % DGTs (original, its adjoint, and one with the differentiated window)
    G = @(x) FDGT(x, tight_win, sigIdx, M, rotIdx, zeroPhaseFlag);
    G_adj = @(u) invFDGT(u, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*w;
    G_diff = @(x) FDGT(x, diff_win, sigIdx, M, rotIdx, zeroPhaseFlag);
    
    % function to calculate the instantaneous frequency of the input signal
    omega = @(x) calcInstFreq(G(x), G_diff(x), M, w, rotateFlag);
    
    % operator to correct phase rotation and its adjoint
    R = @(z, omega) instPhaseCorrection(z, omega, a, M);
    R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);
    
    % time-directional difference
    D = @(z) z(:,1:end-1) - z(:,2:end);
    D_adj = @(z) [z(:,1), (z(:,2:end) - z(:,1:end-1)), -z(:,end)];
    
    % iPC-DGT
    hatG = @(x, omega) D(R(G(x), omega));
    hatG_adj = @(u, omega) G_adj(R_adj(D_adj(u), omega));
    
    
    %%
    
    soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);
    param.proj = @(x) projGamma(x, mask, insig.*mask);
    
    if strcmp(param.type,'B')
    
        sigma = paramsolver.sigma;
        lambda = paramsolver.lambda;
        param.prox = @(z) soft(z, lambda/sigma);
    
        omega_y = omega(insig);
        param.L = @(x) hatG(x, omega_y);
        param.L_adj = @(u) hatG_adj(u, omega_y);
    
        paramsolver.x0 = insig;
        paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
    
        [outsig, snr_procedure] = CP(param, paramsolver, oracle, mask);
    
    
    elseif strcmp(param.type,'Bora')
    
        sigma = paramsolver.sigma;
        lambda = paramsolver.lambda;
        param.prox = @(z) soft(z, lambda/sigma);
    
        omega_y_tilde = omega(oracle);
        param.L = @(x) hatG(x, omega_y_tilde);
        param.L_adj = @(u) hatG_adj(u, omega_y_tilde);
    
        paramsolver.x0 = insig;
        paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
    
        [outsig, snr_procedure] = CP(param, paramsolver, oracle, mask);
    
    
    elseif strcmp(param.type,'R')
    
        wts = 1;
        sigma = paramsolver.sigma;
        lambda = paramsolver.lambda;
    
        omega_y = omega(insig);
        param.L = @(x) hatG(x, omega_y);
        param.L_adj = @(u) hatG_adj(u, omega_y);
    
        paramsolver.x0 = insig;
        paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
        x_old = insig;
    
        snr_procedure = NaN(paramsolver.I, paramsolver.J);
    
        for j = 1:paramsolver.J
    
            param.prox = @(z) soft(z, wts.*lambda/sigma);
    
            [x_hat, snr_procedure(:, j)] = CP(param, paramsolver, oracle, mask);
            
            if norm(x_old - x_hat) < paramsolver.epsilon
                break
            end
    
            denom = movmean(abs(G(x_hat)),[0,1],2);
            wts = 1./(denom + 1e-3);
            wts = wts(:, 1:end-1);
    
            x_old = x_hat;
    
        end
    
        outsig = x_hat;
    
    
    elseif strcmp(param.type,'Rora')
        
        wts = 1;
        sigma = paramsolver.sigma;
        lambda = paramsolver.lambda;
    
        omega_y_tilde = omega(oracle);
        param.L = @(x) hatG(x, omega_y_tilde);
        param.L_adj = @(u) hatG_adj(u, omega_y_tilde);
    
        paramsolver.x0 = insig;
        paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
        x_old = insig;
    
        snr_procedure = NaN(paramsolver.I, paramsolver.J);
    
        for j = 1:paramsolver.J
    
            param.prox = @(z) soft(z, wts.*lambda/sigma);
    
            [x_hat, snr_procedure(:, j)] = CP(param, paramsolver, oracle, mask);
            
            if norm(x_old - x_hat) < paramsolver.epsilon
                break
            end
    
            denom = movmean(abs(G(x_hat)),[0,1],2);
            wts = 1./(denom + 1e-3);
            wts = wts(:, 1:end-1);
    
            x_old = x_hat;
    
        end
    
        outsig = x_hat;
    
    
    elseif strcmp(param.type,'UR')
    
        wts = 1;
        sigma = paramsolver.sigma;
        lambda = paramsolver.lambda;
    
        omega_y = omega(insig);
        param.L = @(x) hatG(x, omega_y);
        param.L_adj = @(u) hatG_adj(u, omega_y);
    
        paramsolver.x0 = insig;
        paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
        x_old = insig;
    
        snr_procedure = NaN(paramsolver.I, paramsolver.J);
    
        for j = 1:paramsolver.J
    
            param.prox = @(z) soft(z, wts.*lambda/sigma);
    
            [x_hat, snr_procedure(:, j)] = CP(param, paramsolver, oracle, mask);
            
            if norm(x_old - x_hat) < paramsolver.epsilon
                break
            end
    
            denom = movmean(abs(G(x_hat)),[0,1],2);
            wts = 1./(denom + 1e-3);
            wts = wts(:, 1:end-1);
    
            omega_x_hat = omega(x_hat);
            param.L = @(x) hatG(x, omega_x_hat);
            param.L_adj = @(u) hatG_adj(u, omega_x_hat);
    
            x_old = x_hat;
    
        end
    
        outsig = x_hat;
    
    
    elseif strcmp(param.type,'U')
    
        sigma = paramsolver.sigma;
        lambda = paramsolver.lambda;
        param.prox = @(z) soft(z, lambda/sigma);
    
        omega_y = omega(insig);
        param.L = @(x) hatG(x, omega_y);
        param.L_adj = @(u) hatG_adj(u, omega_y);
    
        paramsolver.x0 = insig;
        paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));
        x_old = insig;
    
        snr_procedure = NaN(paramsolver.I, paramsolver.J);
    
        for j = 1:paramsolver.J
    
            [x_hat, snr_procedure(:, j)] = CP(param, paramsolver, oracle, mask);
            
            if norm(x_old - x_hat) < paramsolver.epsilon
                break
            end
    
            omega_x_hat = omega(x_hat);
            param.L = @(x) hatG(x, omega_x_hat);
            param.L_adj = @(u) hatG_adj(u, omega_x_hat);
    
            x_old = x_hat;
    
        end
    
        outsig = x_hat;
    
        
    end
end
%% functions 

function x = projGamma(x,mask , corrupted)
    x(mask) = corrupted(mask);
end