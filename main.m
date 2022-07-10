% N_R = 128;
N_R_list = [128, 256, 512, 1024, 2048, 4096];
% seed_value = 1;
seed_value_list = [1, 37, 53];
K = 2;
% Rth = 0*ones(K, 1);
Rth_list = [0, 2].*ones(K, 1);
totalJointIters = 20;
pp = 0;
for seed_value = seed_value_list
    fprintf('seed value %d\n', seed_value);
    for N_R = N_R_list
        fprintf('Number of REs %d\n', N_R);
        for Rth = Rth_list
            fprintf('Rth %d\n', Rth(1));
            funJointOptimisation_SCAref_WMMSE(seed_value, N_R, K, Rth, totalJointIters, pp)
            funJointOptimisation_WMMSE_SCAref(seed_value, N_R, K, Rth, totalJointIters, pp)
            funJointOptimisation_SCAred_WMMSE(seed_value, N_R, K, Rth, totalJointIters, pp)
            funJointOptimisation_WMMSE_SCAred(seed_value, N_R, K, Rth, totalJointIters, pp)
            funJointOptimisation_SCA_WMMSE(seed_value, N_R,    K, Rth, totalJointIters, pp)
        end
    end
end
