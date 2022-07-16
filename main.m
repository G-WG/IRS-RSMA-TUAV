% N_R = 128;
N_R_list = [128, 512, 2048];
% seed_value = 1;
seed_value_list = [110:120];
K = 2;
% Rth = 0*ones(K, 1);
Rth_list = [0, 8].*ones(K, 1);
totalJointIters = 100;
precoderIC = [1, 2]; % {'AMBF-RZF', 'SVD-MRT'} 
pp = 0;
for seed_value = seed_value_list
    fprintf('seed value %d\n', seed_value);
    for prec=precoderIC
        fprintf('precoderIC %d\n', prec);
        for N_R = N_R_list
            fprintf('Number of REs %d\n', N_R);
            for Rth = Rth_list
                if prec == 2
                    Rth = Rth/4;
                end
                fprintf('Rth %d\n', Rth(1));
    %             funJointOptimisation_SCAref_WMMSE(seed_value, N_R, K, Rth, totalJointIters, pp)
                funJointOptimisation_WMMSE_SCAref(seed_value, N_R, K, Rth, totalJointIters, pp, prec)
    %             funJointOptimisation_SCAred_WMMSE(seed_value, N_R, K, Rth, totalJointIters, pp)
                funJointOptimisation_WMMSE_SCAred(seed_value, N_R, K, Rth, totalJointIters, pp, prec)
    %             funJointOptimisation_SCA_WMMSE(seed_value, N_R,    K, Rth, totalJointIters, pp)
            end
        end
    end
end
