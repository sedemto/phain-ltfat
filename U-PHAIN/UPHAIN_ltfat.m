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

    end
    outsig = x_hat;
end
