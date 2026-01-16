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
% parts of this function are copied from PHAINmain.m in [1]

function [outsig, snr_procedure] = UPHAIN_ltfat(insig, mask, param, paramsolver, oracle)

% param
%   .a
%   .M
%   .w
%   .type 


% paramsolver
%   .sigma .... step size (can be ommited when using soft thresholding)
%   .tau ...... step size
%   .alpha .... relaxation paramter
%   .lambda ... threshold
%   .epsilon .. stop criterion
%   .x0 ....... initial value of primal variable
%   .u0 ....... initial value of dual variable
%   .I ........ number of inner iterations
%   .J ........ number of outer iterations


%%
% note that when using soft sigma can be ommited
lambda = paramsolver.lambda;
soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);

% define proximal operators
param.proj = @(x) x.*(1-mask) + insig.*mask;
param.prox = @(z) soft(z, lambda);

if strcmp(param.type,'U')
    % calculate instantaneous frequency 
    omega_y = param.omega(insig);
    
    % setup operator R and its adjoint R*
    phaseCor = exp(-1i*param.phaseCor(omega_y));
    invPhaseCor = exp(1i*param.phaseCor(omega_y));
    param.R = @(z) phaseCor.*z;
    param.R_adj =  @(z) invPhaseCor.*z;

    % combine operators inside the L1 norm into one: L=DRG, L*=G*R*D*
    param.L = @(x) param.D(param.R(param.G(x)));
    param.L_adj = @(u) param.G_adj(param.R_adj(param.D_adj(u)));
    
    % initialize parameters x0 and u0
    paramsolver.x0 = insig;
    paramsolver.u0 = zeros(size(param.L(zeros(length(insig), 1))));

    x_old = insig;

    snr_procedure = NaN(paramsolver.I, paramsolver.J);
    for j = 1:paramsolver.J

        [x_hat, snr_procedure(:, j)] = CP(param, paramsolver, oracle, mask);
        
        if norm(x_old - x_hat) < paramsolver.epsilon
            break
        end
        
        % caluculate new instFreq 
        omega_x_hat = param.omega(x_hat);
        
        % redefine R, R*, L and L* with new instFreq
        phaseCor = exp(-1i*param.phaseCor(omega_x_hat));
        invPhaseCor = exp(1i*param.phaseCor(omega_x_hat));
        param.R = @(z) phaseCor.*z;
        param.R_adj =  @(z) invPhaseCor.*z;

        param.L = @(x) param.D(param.R(param.G(x)));
        param.L_adj = @(u) param.G_adj(param.R_adj(param.D_adj(u)));
        
        x_old = x_hat;
        
        % set previous solution as the new initialization for CP
        if param.updateInputCP
            paramsolver.x0 = x_hat;
        end
    end
    
    outsig = x_hat;
end