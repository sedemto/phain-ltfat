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
% this is a copy of the Matlab function CP.m from [1]

function [x, snr_procedure] = CP(param, paramsolver, oracle, mask)
%% initialization

x = paramsolver.x0;
u = paramsolver.u0;

tau = paramsolver.tau;
alpha = paramsolver.alpha;

snr_procedure = NaN(paramsolver.I, 1);
%% iteration

for i = 1:paramsolver.I

    p = param.proj(x - tau*param.L_adj(u));
    v = u + param.L(2*p - x);
    q = v - param.prox(v);
    
    x = x + alpha*(p - x);
    u = u + alpha*(q - u);
    
    snr_procedure(i) = snr(oracle(~mask),oracle(~mask)- x(~mask));
end
