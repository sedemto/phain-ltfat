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

hatG = @(x, omega) param.D(param.R(param.G(x), omega));
hatG_adj = @(u, omega) param.G_adj(param.R_adj(param.D_adj(u), omega));

%%

soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);
param.proj = @(x) x.*(1-mask) + insig.*mask;

if strcmp(param.type,'U')
    sigma = paramsolver.sigma;
    lambda = paramsolver.lambda;
    param.prox = @(z) soft(z, lambda/sigma);

    omega_y = param.omega(insig);
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

        omega_x_hat = param.omega(x_hat);
        param.L = @(x) hatG(x, omega_x_hat);
        param.L_adj = @(u) hatG_adj(u, omega_x_hat);
        
        x_old = x_hat;
        if param.updateInputCP
            paramsolver.x0 = x_hat;
        end
    end
    
    outsig = x_hat;
end



