N_R = 128;
seed_value = 1;
K = 2;
Rth = 0*ones(K, 1);
totalJointIters = 20;
pp = 0;

funJointOptimisation_SCAref_WMMSE(seed_value, N_R, K, Rth, totalJointIters, pp)
funJointOptimisation_WMMSE_SCAref(seed_value, N_R, K, Rth, totalJointIters, pp)
funJointOptimisation_SCAred_WMMSE(seed_value, N_R, K, Rth, totalJointIters, pp)
funJointOptimisation_WMMSE_SCAred(seed_value, N_R, K, Rth, totalJointIters, pp)
funJointOptimisation_SCA_WMMSE(seed_value, N_R,    K, Rth, totalJointIters, pp)
