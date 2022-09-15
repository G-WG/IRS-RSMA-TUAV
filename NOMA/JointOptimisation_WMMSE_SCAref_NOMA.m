%%
% JointOptimisation_WMMSE_SCAref_NOMA
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
% - 4.0
%   - Modifications for NOMA scheme. It based on encoding W1 to s1 and W2
%   to s12, and no power to s2 (or vice-versa).


%%

tau = 0.8; h_B = 30;
q_B = [0, 0, h_B]';
initialise_params_differentScenarios2nd_v4_NOMA;
% weights
u_k = ones(K, 1);
% Threshold rate
% Rth = 0*ones(K, 1);
epsilonNOMA = 1e-6;
epsilonPS = 1e-6;
maxIter = 1;

NOMAoptCell = {};
PSoptCell = {};

NOMAtimeElapsed = 0;
SCAtimeElapsed = 0;
if precoderIC == 1
    name = sprintf('NOMA_JointOptimisation_NOMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_AMBF_RZF_TUAV_%d_%d_%d', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R, q_TUAV(1), q_TUAV(2), q_TUAV(3));
else
    name = sprintf('NOMA_JointOptimisation_NOMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_SVD_MRT_TUAV_%d_%d_%d', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R, q_TUAV(1), q_TUAV(2), q_TUAV(3));
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

while loop

%% NOMA MMSE
if fileID
    fprintf(fileID, '\n');
    fprintf(fileID, '##########  NOMA  ##########\n');
    fprintf(fileID, '\n');
end
Theta = diag(s);
h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';

try
    [NOMAoptStructureWMMSE, NOMAloopWMMSE, NOMAiterationsWMMSE, timeElapsedWMMSE] = ...
            Optimisation_Update_NOMAparameters_WMMSE(h_ov_k, K, ...
            varianceNoise, p_c_IC, p_k_IC, Pt, N_T, Rth, u_k, common_rates, ...
            epsilonNOMA, maxIter, fileID, UE_NOMA_common, UE_NOMA_private);
    p_c_IC = NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Pc;
    p_k_IC = NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Pk;
    NOMAoptCell{iter} = {NOMAoptStructureWMMSE(1:NOMAiterationsWMMSE), NOMAloopWMMSE, NOMAiterationsWMMSE, timeElapsedWMMSE};
    common_rates = NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Cc;
    NOMAtimeElapsed = NOMAtimeElapsed + timeElapsedWMMSE;
    if wsr < NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Rov
        wsr = NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Rov;
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
        NOMAconvergence =  abs(NOMAoptCell{iter}{1}(end).Rov - NOMAoptCell{iter-1}{1}(end).Rov) <= epsilonNOMA;
        PSconvergence   =  abs(PSoptCell{iter}{1}(end).Rov - PSoptCell{iter-1}{1}(end).Rov) <= epsilonPS;
        if NOMAconvergence && PSconvergence
            loop = 0;
        end
    
    end
    if fileID
        fprintf('iter: %d, PS-Rov=%.8f, NOMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Rov)
        fprintf(fileID, 'iter: %d, PS-Rov=%.8f, NOMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Rov);
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
    fprintf(fileID, 'NOMA time Elapsed: %.8f [s], SCA time Elapsed: %.8f [s], total Time: %.8f [s]\n', NOMAtimeElapsed, SCAtimeElapsed, NOMAtimeElapsed+SCAtimeElapsed);
    fprintf(fileID, '\n');
    fprintf(fileID, '##########  END  ##########\n');
    fprintf(fileID, '\n');
    fclose(fileID);
end

if matfileHandler
    name_matfile = sprintf('logs/JointOptimisation/log_seed_%d_%s.mat', seed_value, name);
    save(name_matfile)
end