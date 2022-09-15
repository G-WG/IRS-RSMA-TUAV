function [optStructure, loop, jj_phaseshift, timeElapsed] = ...
        Optimisation_Update_PhaseShiftref_SDMA(s, h_T_U_PL, h_R_U_PL, G, K, ...
        varianceNoise, p_c_IC, p_k_IC, N_R, Rth, u_k, common_rates, epsilon, maxIter, fileID)


% [optStructure, loop, jj_phaseshift] = ...
%         Optimisation_Update_PhaseShiftref_SDMA(s, h_T_U_PL, h_R_U_PL, G, K, ...
%         varianceNoise, p_c_IC, p_k_IC, N_R, Rth, u_k, epsilon, maxIter, fileID)
%
% This functions performs iterations over the phase shift optimisation.
% 
% Input parameters:
% - s : Initial condition for the phase shift. Size of [N_R, 1]
% - h_T_U_PL : Channel from the TUAV to the user.
% - h_R_U_PL : Channel from the IRS to the user.
% - G : Channel from the TUAV to the IRS.
% - K :  number of users.
% - varianceNoise : Noise variance in Watts.
% - p_c_IC : Precoder vector for the common stream.
% - p_k_IC : Precoder matrix for the private stream. Size of [N_T, K], N_T
% - N_R :  number of reflecting elements.
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
%     optStructure_.Rp = zeros(K, 1);
%     optStructure_.Rc = zeros(K, 1);
%     optStructure_.Rov = zeros(1);
%     optStructure_.phaseshifts = zeros(N_R, 1);
%     optStructure_.PSgt1 = zeros(1);
%     optStructure_.optimalValue = zeros(1);
% - loop: either 1 or 0. If it is 0, then the optimisation converge,
% otherwise, the optimisation reach the limit of iterations.
% - jj_phaseshift: the number of iterations to converge.
%
% A typical entry could be:
% fileID = 0;
% initialise_params;
% u_k = ones(K, 1);
% Rth = 0*ones(K, 1);
% epsilonPS = 1e-4;
% maxIter = 1000;
% - [PSoptStructure, PSloop, PSiterations] = Optimisation_Update_PhaseShiftred(s, h_T_U_PL, ...
%         h_R_U_PL, G, K, varianceNoise, p_c_IC, p_k_IC, N_R, Rth, u_k, ...
%         epsilonPS, maxIter, fileID);
%
% Nested Functions:
% - compute_rates.m
% - CVX_opt_phase_shift_ref_SDMA.m
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/06/30

%% Structure
optStructure_.iter = 0;
optStructure_.Rp = zeros(K, 1);
% optStructure_.Rc = zeros(K, 1);
optStructure_.Rov = zeros(1);
optStructure_.phaseshifts = zeros(N_R, 1);
optStructure_.PSgt1 = zeros(1);
optStructure_.optimalValue = zeros(1);
optStructure = repmat(optStructure_, maxIter+1, 1);
% initial condition for the optimal value
optStructure(1).optimalValue = 1e6;
%% Initialise parameters

Theta = diag(s);
h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
[rate_c, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
% common_rates = min(rate_c).*ones(K, 1)/K;
% Lambda_c = 2.^sum(common_rates)-1;
Lambda_p = 2.^(Rth-common_rates)-1;
s_IC = s;
eta_p_IC = (2.^rate_kp-1);
beta_p_IC = zeros(K, 1);
C = 1e4;
h_ov_k = (h_T_U_PL' + h_R_U_PL' * diag(s) * G)';
for k = 1:K
    interference_term = 0;
    for jj_phaseshift = 1:K
        if k ~= jj_phaseshift
            interference_term = interference_term + ...
            power(abs(h_ov_k(:, k)'*p_k_IC(:, jj_phaseshift)), 2);
        end
    end
    beta_p_IC(k) = interference_term + varianceNoise;
end

loop = 1;

if fileID
    fprintf(fileID, 'Initial Conditions: ');
    for ii = 1:K
        fprintf(fileID, 'R%dp %.8f, ', ii, rate_kp(ii));
    end
%     for ii = 1:K
%         fprintf(fileID, 'R%dc %.8f, ', ii, rate_c(ii));
%     end
    fprintf(fileID, 'Rov = %.8f\n', u_k'*(rate_kp+min(rate_c)/K));
end

% common_rates = min(rate_c)/K;
optStructure(1).Rp = rate_kp;
% optStructure(1).Rc = rate_c;
% optStructure(1).Cc = common_rates;
optStructure(1).Rov = u_k'*(rate_kp);
optStructure(1).P = [trace(p_c_IC*p_c_IC')]';

%% Iteration

% fprintf(fileID, '\n');
% fprintf(fileID, '##########  PS  ##########\n');
% fprintf(fileID, '\n');

jj_phaseshift = 1;
tic;
while loop

    [opt_val, s, eta_p, beta_p, status] = CVX_opt_phase_shift_ref_SDMA(s_IC, ...
       eta_p_IC, beta_p_IC, C, K, N_R, h_T_U_PL, h_R_U_PL, G, p_k_IC, ...
       varianceNoise, Lambda_p);
    
    s_IC = s;
    eta_p_IC = eta_p;
    beta_p_IC = beta_p;
    Theta = diag(s);
    h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
    [rate_c, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
    
    if fileID
        fprintf(fileID, 'opt value: %.8f at %3d, PS gt1 %4d, ', opt_val, jj_phaseshift, sum(round(abs(s), 8) > 1));
        for ii = 1:K
            fprintf(fileID, 'R%dp %.8f, ', ii, rate_kp(ii));
        end
%         for ii = 1:K
%             fprintf(fileID, 'R%dc %.8f, ', ii, rate_c(ii));
%         end
        fprintf(fileID, 'Rov = %.8f, ', u_k'*(rate_kp+min(rate_c)/K));
        fprintf(fileID, 'Status: %s\n', status);
    end    
    

    optStructure(jj_phaseshift+1).optimalValue = opt_val;
    optStructure(jj_phaseshift+1).iter = jj_phaseshift;
    optStructure(jj_phaseshift+1).Rp = rate_kp;
%     optStructure(jj_phaseshift+1).Rc = rate_c;
    optStructure(jj_phaseshift+1).Rov = u_k'*(rate_kp);
    optStructure(jj_phaseshift+1).phaseshifts = s;
    optStructure(jj_phaseshift+1).PSgt1 = sum(abs(s) > 1);

    if abs(optStructure(jj_phaseshift+1).optimalValue - optStructure(jj_phaseshift).optimalValue ) < epsilon
        loop = 0;
    end
    if jj_phaseshift >= maxIter
        jj_phaseshift = jj_phaseshift + 1;
        break
    end
    jj_phaseshift = jj_phaseshift + 1;
end
timeElapsed = toc;
if fileID
    fprintf(fileID, 'TimeElapsed = %.8f\n', timeElapsed);
end



end