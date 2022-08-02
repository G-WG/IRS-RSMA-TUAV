% clear all; 
% seed_value = 1; N_R = 128; totalJointIters = 20;
% pp = 0;
% fileID = 1;
tau = 0;
initialise_params_TUAVv2;
% weights
u_k = ones(K, 1);
% Threshold rate
% Rth = 0*ones(K, 1);
epsilonNOMA = 1e-6;
epsilonPS = 1e-4;
maxIter = 1000;

NOMAoptCell = {};
PSoptCell = {};

NOMAtimeElapsed = 0;
SCAtimeElapsed = 0;
if precoderIC == 1
    name = sprintf('JointOptimisation_NOMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_AMBF_RZF', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R);
else
    name = sprintf('JointOptimisation_NOMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_SVD_MRT', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R);
end
fileID = fopen(sprintf('logs/JointOptimisation/log_seed_%d_%s.txt', seed_value, name),'w');
iter = 1;
common_rates = min(rate_c).*ones(K, 1)/K;
loop = 1;
while loop

%% NOMA MMSE
fprintf(fileID, '\n');
fprintf(fileID, '##########  NOMA  ##########\n');
fprintf(fileID, '\n');
Theta = diag(s);
h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';

try
    [NOMAoptStructureWMMSE, NOMAloopWMMSE, NOMAiterationsWMMSE, timeElapsedWMMSE] = ...
            Optimisation_Update_NOMAparameters_WMMSE(h_ov_k, K, ...
            varianceNoise, p_c_IC, Pt, N_T, Rth, u_k, epsilonNOMA, maxIter, fileID);

    p_c_IC = NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Pc;
    p_k_IC = zeros(N_T, K);
    NOMAoptCell{iter} = {NOMAoptStructureWMMSE(1:NOMAiterationsWMMSE), NOMAloopWMMSE, NOMAiterationsWMMSE, timeElapsedWMMSE};
    common_rates = NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Cc;
    NOMAtimeElapsed = NOMAtimeElapsed + timeElapsedWMMSE;
catch ME
    disp(ME)
    fprintf(fileID, 'Infeasible.\n');
end
%% Phase shift

fprintf(fileID, '\n');
fprintf(fileID, '##########  PS  ##########\n');
fprintf(fileID, '\n');

try
    [PSoptStructureSCAred, PSloopSCAred, PSiterationsSCAred, PStimeElapsedSCAred] = ...
            Optimisation_Update_PhaseShiftref(s, h_T_U_PL, h_R_U_PL, G, K, ...
        varianceNoise, p_c_IC, p_k_IC, N_R, Rth, u_k, common_rates, epsilonPS, maxIter, fileID);
    s = PSoptStructureSCAred(PSiterationsSCAred).phaseshifts;
    PSoptCell{iter} = {PSoptStructureSCAred(1:PSiterationsSCAred), PSloopSCAred, PSiterationsSCAred, PStimeElapsedSCAred};
    SCAtimeElapsed = SCAtimeElapsed + PStimeElapsedSCAred;
catch ME
    disp(ME)
    fprintf(fileID, 'Infeasible.\n');
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
    fprintf('iter: %d, PS-Rov=%.8f, NOMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Rov);
    fprintf(fileID, 'iter: %d, PS-Rov=%.8f, NOMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, NOMAoptStructureWMMSE(NOMAiterationsWMMSE).Rov);
catch ME
    disp(ME)
    fprintf(fileID, 'CVX non-return.\n');
    break;
end


if iter > totalJointIters
    loop = 0;
end


iter = iter + 1;



%% save workspace and close log
end
fprintf('NOMA time Elapsed: %.8f [s], SCA time Elapsed: %.8f [s], total Time: %.8f [s]\n', NOMAtimeElapsed, SCAtimeElapsed, NOMAtimeElapsed+SCAtimeElapsed)
fprintf(fileID, '\n');
fprintf(fileID, '##########  END  ##########\n');
fprintf(fileID, '\n');
fclose(fileID);
name_matfile = sprintf('logs/JointOptimisation/log_seed_%d_%s.mat', seed_value, name);
save(name_matfile)