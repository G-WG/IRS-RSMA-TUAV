%%
% JointOptimisation_SDMA_WMMSE_SCAref.m
% author: Maximiliano Rivera
% date: 09-02-2022
% version: 2.0.0
%   - channel gain: initialise_params_TUAVv2
% version: 2.0.1 
%   - Using: initialise_params_differentScenarios2nd_v4.m 

%%

tau = 1; % 1 for only private precoder
h_B = 30;
q_B = [0, 0, h_B]';

% initialise_params_TUAVv2;
initialise_params_differentScenarios2nd_v4;
% weights
u_k = ones(K, 1);
% Threshold rate
% Rth = 0*ones(K, 1);
epsilonSDMA = 1e-6;
epsilonPS = 1e-6;
maxIter = 1;

SDMAoptCell = {};
PSoptCell = {};

SDMAtimeElapsed = 0;
SCAtimeElapsed = 0;
if precoderIC == 1
    name = sprintf('JointOptimisation_SDMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_AMBF_RZF', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R);
else
    name = sprintf('JointOptimisation_SDMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_SVD_MRT', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R);
end
fileID = fopen(sprintf('logs/JointOptimisation/log_seed_%d_%s.txt', seed_value, name),'w');
iter = 1;
common_rates = min(rate_c).*ones(K, 1)/K;
loop = 1;
if fileID
    fprintf(fileID, 'q_TUAV = [%.0f, %.0f, %.0f]\n', q_TUAV(1), q_TUAV(2), q_TUAV(3));
end

while loop

%% SDMA MMSE
fprintf(fileID, '\n');
fprintf(fileID, '##########  SDMA  ##########\n');
fprintf(fileID, '\n');
Theta = diag(s);
h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';

try

    [SDMAoptStructureWMMSE, SDMAloopWMMSE, SDMAiterationsWMMSE, timeElapsedWMMSE] = ...
        Optimisation_Update_SDMAparameters_WMMSE(h_ov_k, K, ...
        varianceNoise, p_k_IC, Pt, N_T, Rth, u_k, epsilonSDMA, maxIter, fileID);

    p_k_IC = SDMAoptStructureWMMSE(SDMAiterationsWMMSE).Pk;
    SDMAoptCell{iter} = {SDMAoptStructureWMMSE(1:SDMAiterationsWMMSE), SDMAloopWMMSE, SDMAiterationsWMMSE, timeElapsedWMMSE};
    SDMAtimeElapsed = SDMAtimeElapsed + timeElapsedWMMSE;
    if wsr < SDMAoptStructureWMMSE(SDMAiterationsWMMSE).Rov
        wsr = SDMAoptStructureWMMSE(SDMAiterationsWMMSE).Rov;
    end
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
            Optimisation_Update_PhaseShiftref_SDMA(s, h_T_U_PL, h_R_U_PL, G, K, ...
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
        SDMAconvergence =  abs(SDMAoptCell{iter}{1}(end).Rov - SDMAoptCell{iter-1}{1}(end).Rov) <= epsilonSDMA;
        PSconvergence   =  abs(PSoptCell{iter}{1}(end).Rov - PSoptCell{iter-1}{1}(end).Rov) <= epsilonPS;
        if SDMAconvergence && PSconvergence
            loop = 0;
        end
    
    end
    fprintf('iter: %d, PS-Rov=%.8f, SDMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, SDMAoptStructureWMMSE(SDMAiterationsWMMSE).Rov);
    fprintf(fileID, 'iter: %d, PS-Rov=%.8f, SDMA-Rov=%.8f\n', iter, PSoptStructureSCAred(PSiterationsSCAred).Rov, SDMAoptStructureWMMSE(SDMAiterationsWMMSE).Rov);
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
if fileID
    fprintf(fileID, 'SDMA time Elapsed: %.8f [s], SCA time Elapsed: %.8f [s], total Time: %.8f [s]\n', SDMAtimeElapsed, SCAtimeElapsed, SDMAtimeElapsed+SCAtimeElapsed)
    fprintf(fileID, '\n');
    fprintf(fileID, '##########  END  ##########\n');
    fprintf(fileID, '\n');
    fclose(fileID);
end
if matfileHandler
    name_matfile = sprintf('logs/JointOptimisation/log_seed_%d_%s.mat', seed_value, name);
    save(name_matfile)
end