function wsr = funJointOptimisation_SDMA_WMMSE_SCAref(seed_value, N_R, K, Rth, totalJointIters, pp, precoderIC, q_TUAV)
%%
% funJointOptimisation_SDMA_WMMSE_SCAref(seed_value, N_R, K, Rth, totalJointIters, pp, precoderIC)
% author: Maximiliano Rivera
% date: 09-02-2022
% version: 
% - 1.0
% - 2.0: TUAV placement, filename, matfile handler
disp('Start: JointOptimisation_SDMA_WMMSE_SCAref')
matfileHandler = 1;
fileIDHandler = 1;
JointOptimisation_SDMA_WMMSE_SCAref;
disp('finish JointOptimisation_SDMA_WMMSE_SCAref')
end