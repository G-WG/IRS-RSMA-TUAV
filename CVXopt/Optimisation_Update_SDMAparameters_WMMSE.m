function [optStructure, loop, jj, timeElapsed] = ...
        Optimisation_Update_SDMAparameters_WMMSE(h_ov_k, K, ...
        varianceNoise, p_k_IC, Pt, N_T, Rth, u_k, epsilon, maxIter, fileID)

% [optStructure, loop, jj] = ...
%         Optimisation_Update_SDMAparameters_WMMSE(h_ov_k, K, ...
%         varianceNoise, p_k_IC, Pt, N_T, Rth, u_k, epsilon, maxIter, fileID)
%
% This functions performs iterations over the RSMA parameters optimisation 
% using SCA.
% 
% Input parameters:
% - h_ov_k : Overall Channel. Size of [N_T, K]
% - K :  number of users.
% - varianceNoise : Noise variance in Watts.
% - p_k_IC : Precoder matrix for the private stream. Size of [N_T, K], N_T
% - Pt : Maximum transmit power.
% - N_T :  Number of transmit antennas.
% - Rth : rate threshold. size of [K, 1]
% - u_k : weights for the WSR. size of [K, 1]
% - epsilon : convergence resolution.
% - maxIter : Maximum number of iterations.
% - fileID : file variable from fopen function. If file=1, then the log is
% printed. If file=0, there is no log of the file displayed.
%
% Output parameters:
% - optStructure: Structure containing all the relevan information about
% the optimisation. It is built by:
%     optStructure_.iter = 0;
%     optStructure_.Rp = zeros(K, 1); -- Private rate
%     optStructure_.Rc = zeros(K, 1); -- Common rate
%     optStructure_.Rov = zeros(1);  -- Overall rate
%     optStructure_.Cc = zeros(K, 1); -- Common allocated rate
%     optStructure_.P = zeros(2, 1); -- Power used in both common and
%                                       private precoders.
%     optStructure_.optimalValue = zeros(1);
% - loop: either 1 or 0. If it is 0, then the optimisation converge,
% otherwise, the optimisation reach the limit of iterations.
% - jj: the number of iterations to converge.
%
% A typical entry could be:
% fileID = 0;
% initialise_params;
% u_k = ones(K, 1);
% Rth = 0*ones(K, 1);
% epsilonPS = 1e-4;
% maxIter = 1000;
% - [optStructure, loop, jj] = ...
%         Optimisation_Update_RSMAparameters_SCA(h_ov_k, K, ...
%         varianceNoise, p_c_IC, p_k_IC, Pt, N_T, Rth, u_k, epsilon, maxIter, fileID)
%
% Nested Functions:
% - compute_rates.m
% - CVX_opt_SDMA_parameters_WMMSE.m
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/06/30

%% Structure
optStructure_.iter = 0;
optStructure_.Rp = zeros(K, 1);
optStructure_.Rov = zeros(1);
optStructure_.P = zeros(1, 1); % (||Pk||)
optStructure_.Pk = zeros(N_T, K);
optStructure_.optimalValue = zeros(1);
optStructure = repmat(optStructure_, maxIter, 1);
% initial condition for the optimal value
optStructure(1).optimalValue = 1e6;
%% Initialise parameters
loop = 1;
p_c_IC = zeros(N_T, 1);
[~, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
if fileID
    fprintf(fileID, 'Initial Conditions: ');
    for ii = 1:K
        fprintf(fileID, 'R%dp %.8f, ', ii, rate_kp(ii));
    end
    fprintf(fileID, 'Rov = %.8f, ', u_k'*(rate_kp));
    fprintf(fileID, '||Pk|| %.8f\n', trace(p_k_IC*p_k_IC'));
end


optStructure(1).Rp = rate_kp;
optStructure(1).Rov = u_k'*(rate_kp);
optStructure(1).P = [trace(p_k_IC*p_k_IC')]';
optStructure(1).Pk = p_k_IC;

%% CVX RSMA param


jj = 1;
tic;
while loop

    [opt_val, p_k] = ...
        CVX_opt_SDMA_parameters_WMMSE(p_k_IC,...
        K, N_T, u_k, Rth, Pt, h_ov_k, varianceNoise);
    

    p_k_IC = p_k;
    [~, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k);

    if fileID
        fprintf(fileID, 'opt value: %.8f at %3d, ', opt_val, jj);
        for ii = 1:K
            fprintf(fileID, 'R%dp %.8f, ', ii, rate_kp(ii));
        end
        fprintf(fileID, 'Rov = %.8f, ', u_k'*(rate_kp));
        fprintf(fileID, '||Pk|| %.8f\n',trace(p_k*p_k'));
    end
   

    optStructure(jj+1).optimalValue = opt_val;
    optStructure(jj+1).iter = jj;
    optStructure(jj+1).Rp = rate_kp;
    optStructure(jj+1).Rov = u_k'*(rate_kp);
    optStructure(jj+1).P = [trace(p_k*p_k')]';
    optStructure(jj+1).Pk = p_k_IC;

    if abs(optStructure(jj+1).optimalValue - optStructure(jj).optimalValue ) < epsilon
        loop = 0;
    end
    if jj >= maxIter
        break
    end
    jj = jj + 1;
end
timeElapsed = toc;
if fileID
    fprintf(fileID, 'TimeElapsed = %.8f\n', timeElapsed);
end

end