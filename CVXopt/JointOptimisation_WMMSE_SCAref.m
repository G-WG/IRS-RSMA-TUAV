%%
% JointOptimisation_WMMSE_SCAref
% author: Maximiliano Rivera
% date: 09-02-2022
% version: 
% - 2.0
%   - Channel gain fixed: initialise_params_TUAVv2
% - 3.0
%   - TUAV channel environment improved: LOS is evaluated through building
%   blockage.
%   - Using: initialise_params_differentScenarios2nd_v4.m
% - 3.1
%   - Add TUAV position, file name and matfile handlers


%%

% clear all; 
% name = 'PSred_RSMA_SCA';
% seed_value = 1; N_R = 128; totalJointIters = 20;
% pp = 0;
% fileID = fopen(sprintf('logs/RSMA_comparison_logs/log_seed_%d_%s.txt', seed_value, name),'w');
% fileID = 1;
tau = 0.8; h_B = 30;
q_B = [0, 0, h_B]';
% q_TUAV = [0, 0, 100]' + q_B;
% initialise_params_TUAVv2;
initialise_params_differentScenarios2nd_v4;
% weights
u_k = ones(K, 1);
% Threshold rate
% Rth = 0*ones(K, 1);
epsilonRSMA = 1e-6;
epsilonPS = 1e-6;
maxIter = 1;

RSMAoptCell = {};
PSoptCell = {};

RSMAtimeElapsed = 0;
SCAtimeElapsed = 0;
if precoderIC == 1
    name = sprintf('JointOptimisation_RSMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_AMBF_RZF_TUAV_%d_%d_%d', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R, q_TUAV(1), q_TUAV(2), q_TUAV(3));
else
    name = sprintf('JointOptimisation_RSMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_SVD_MRT_TUAV_%d_%d_%d', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R, q_TUAV(1), q_TUAV(2), q_TUAV(3));
end
if fileIDHandler
    fileID = fopen(sprintf('logs/JointOptimisation/log_seed_%d_%s.txt', seed_value, name),'w');
else
    fileID = 0;
end
iter = 1;
common_rates = min(rate_c).*ones(K, 1)/K;
loop = 1;
if fileID
    fprintf(fileID, 'q_TUAV = [%.0f, %.0f, %.0f]\n', q_TUAV(1), q_TUAV(2), q_TUAV(3));
end
wsr = 0;
while loop

%% RSMA MMSE
if fileID
    fprintf(fileID, '\n');
    fprintf(fileID, '##########  RSMA  ##########\n');
    fprintf(fileID, '\n');
end
Theta = diag(s);
h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';

try
    [RSMAoptStructureWMMSE, RSMAloopWMMSE, RSMAiterationsWMMSE, timeElapsedWMMSE] = ...
            Optimisation_Update_RSMAparameters_WMMSE(h_ov_k, K, ...
            varianceNoise, p_c_IC, p_k_IC, Pt, N_T, Rth, u_k, common_rates, epsilonRSMA, maxIter, fileID);
    p_c_IC = RSMAoptStructureWMMSE(RSMAiterationsWMMSE).Pc;
    p_k_IC = RSMAoptStructureWMMSE(RSMAiterationsWMMSE).Pk;
    RSMAoptCell{iter} = {RSMAoptStructureWMMSE(1:RSMAiterationsWMMSE), RSMAloopWMMSE, RSMAiterationsWMMSE, timeElapsedWMMSE};
    common_rates = RSMAoptStructureWMMSE(RSMAiterationsWMMSE).Cc;
    RSMAtimeElapsed = RSMAtimeElapsed + timeElapsedWMMSE;
    if wsr < RSMAoptStructureWMMSE(RSMAiterationsWMMSE).Rov
        wsr = RSMAoptStructureWMMSE(RSMAiterationsWMMSE).Rov;
    end
catch ME
    disp(ME)
    if fileID
        fprintf(fileID, 'Infeasible.\n');
    end
end
%% Phase shift
if fileID
    fprintf(fileID, '\n');
    fprintf(fileID, '##########  PS  ##########\n');
    fprintf(fileID, '\n');
end

try
    [PSoptStructureSCAred, PSloopSCAred, PSiterationsSCAred, PStimeElapsedSCAred] = ...
            Optimisation_Update_PhaseShiftref(s, h_T_U_PL, h_R_U_PL, G, K, ...
        varianceNoise, p_c_IC, p_k_IC, N_R, Rth, u_k, common_rates, epsilonPS, maxIter, fileID);
    s = PSoptStructureSCAred(PSiterationsSCAred).phaseshifts;
    PSoptCell{iter} = {PSoptStructureSCAred(1:PSiterationsSCAred), PSloopSCAred, PSiterationsSCAred, PStimeElapsedSCAred};
    SCAtimeElapsed = SCAtimeElapsed + PStimeElapsedSCAred;
    if wsr < PSoptStructureSCAred(PSiterationsSCAred).Rov
        wsr = PSoptStructureSCAred(PSiterationsSCAred).Rov;
    end
catch ME
    disp(ME)
    if fileID
        fprintf(fileID, 'Infeasible.\n');
    end
end

% if PSoptStructureSCAred(PSiterationsSCAred).PSgt1 > 0
%     break
% end




%%
try
    if iter > 1
        RSMAconvergence =  abs(RSMAoptCell{iter}{1}(end).Rov - RSMAoptCell{iter-1}{1}(end).Rov) <= epsilonRSMA;
        PSconvergence   =  abs(PSoptCell{iter}{1}(end).Rov - PSoptCell{iter-1}{1}(end).Rov) <= epsilonPS;
        if RSMAconvergence && PSconvergence
            loop = 0;
        end
    
    end
    if fileID
        fprintf('iter: %d, PS-Rov=%.8f, RSMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, RSMAoptStructureWMMSE(RSMAiterationsWMMSE).Rov)
        fprintf(fileID, 'iter: %d, PS-Rov=%.8f, RSMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, RSMAoptStructureWMMSE(RSMAiterationsWMMSE).Rov);
    end
catch ME
    disp(ME)
    if fileID
        fprintf(fileID, 'CVX non-return.\n');
    end
    break;
end


if iter > totalJointIters
    loop = 0;
end


iter = iter + 1;



%% save workspace and close log
end
if fileID
    fprintf(fileID, 'RSMA time Elapsed: %.8f [s], SCA time Elapsed: %.8f [s], total Time: %.8f [s]\n', RSMAtimeElapsed, SCAtimeElapsed, RSMAtimeElapsed+SCAtimeElapsed);
    fprintf(fileID, '\n');
    fprintf(fileID, '##########  END  ##########\n');
    fprintf(fileID, '\n');
    fclose(fileID);
end

if matfileHandler
    name_matfile = sprintf('logs/JointOptimisation/log_seed_%d_%s.mat', seed_value, name);
    save(name_matfile)
end