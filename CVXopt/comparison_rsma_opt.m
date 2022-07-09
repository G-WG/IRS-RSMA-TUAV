N = 1000; seed_value = 1;initialise_params;
structur4CVX;

CVXstruct_SCA = repmat(CVXstruct, N, 1);
CVXstruct_SCA_ref = repmat(CVXstruct, N, 1);
CVXstruct_WMMSE = repmat(CVXstruct, N, 1);

seed_list = 22;
%%
for seed_value = seed_list
    disp(seed_value)
    initialise_params
    
    [rate_c, rate_kp] = compute_rates(s, h_T_U_PL, h_R_U_PL, G, K, varianceNoise, p_c_IC, p_k_IC);
    % weights
    u_k = ones(K, 1);
    
    % Threshold rate
    Rth = 0*ones(K, 1);
    % Define the common initial rates as the min of both divided by the number
    % of users.
    common_rates = min(rate_c).*ones(K, 1)/K;
    
    Lambda_c = 2.^sum(common_rates)-1;
    Lambda_p = 2.^(Rth-common_rates)-1;
    
    
    s_IC = s;
    eta_p_IC = (2.^rate_kp-1);
    beta_p_IC = zeros(K, 1);
    zeta_p_IC = zeros(K, 1);
    zeta_c_IC = zeros(K, 1);
    C = 1e4;
    for k = 1:K
        interference_term = 0;
        for jj = 1:K
            if k ~= jj
                interference_term = interference_term + ...
                power(abs(h_ov_k(:, k)'*p_k_IC(:, jj)), 2);
            end
        end
    %     k
        beta_p_IC(k) = interference_term + varianceNoise;
        zeta_p_IC(k) = interference_term + varianceNoise;
        zeta_c_IC(k) = interference_term + varianceNoise + pow_abs(h_ov_k(:, k)'*p_k_IC(:, k), 2);
    end
    
    rho_c_IC = 2.^(rate_c)-1; %rho_c_IC = rho_c_IC * 0;
    rho_p_IC = 2.^(rate_kp)-1; %rho_p_IC = rho_p_IC * 0;
    
    epsilon_rsma = 1e-6;
    Theta = diag(s);
    h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
    % fprintf('R1p %.8f, R2p %.8f, ', rate_kp(1), rate_kp(2))
    % fprintf('R1c %.8f, R2c %.8f\n', rate_c(1), rate_c(2))
    %% RSMA SCA
    
    
    
    
    [opt_val, ~, ~, ~, ~, ~, ~, p_k, p_c, common_rates] = ...
            CVX_opt_rsma_parameters(rho_c_IC, rho_p_IC, zeta_c_IC, zeta_p_IC, ...
            p_k_IC, p_c_IC, K, N_T, u_k, Rth, Pt, h_ov_k, varianceNoise);
    
%     fprintf('%4d, %3s: ', seed_value, 'SCA')
%     fprintf('opt value:  %.4f at %3d, ', opt_val, jj);
    [rate_c, rate_kp] = compute_rates(s, h_T_U_PL, h_R_U_PL, G, K, varianceNoise, p_c, p_k);
%     fprintf('R1p %.8f, R2p %.8f, ', rate_kp(1), rate_kp(2))
%     fprintf('R1c %.8f, R2c %.8f, ', rate_c(1), rate_c(2))
%     fprintf('C1c %.8f, C2c %.8f, ', common_rates(1), common_rates(2));
%     fprintf('Rov = %.8f, ', u_k'*(rate_kp+common_rates));
%     fprintf('||Pc|| %.8f, ||Pk|| %.8f\n',trace(p_c*p_c'), trace(p_k*p_k'))
    
    CVXstruct_SCA(seed_value).iter = seed_value;
    CVXstruct_SCA(seed_value).Rp = rate_kp;
    CVXstruct_SCA(seed_value).Rc = rate_c;
    CVXstruct_SCA(seed_value).Cc = common_rates;
    CVXstruct_SCA(seed_value).Rov = u_k'*(rate_kp+common_rates);
    CVXstruct_SCA(seed_value).P = [trace(p_c*p_c'), trace(p_k*p_k')];
    
    %% SCA ref
    
    [opt_val, ~, ~, ~, ~, ~, ~, p_k, p_c, common_rates] = ...
            CVX_opt_rsma_parameters_ref(zeta_c_IC, zeta_p_IC, ...
            p_k_IC, p_c_IC, K, N_T, u_k, Rth, Pt, h_ov_k, varianceNoise);
    
%     fprintf('  SCA ref: ')
%     fprintf('opt value:  %.4f at %3d, ', opt_val, jj);
    [rate_c, rate_kp] = compute_rates(s, h_T_U_PL, h_R_U_PL, G, K, varianceNoise, p_c, p_k);
%     fprintf('R1p %.8f, R2p %.8f, ', rate_kp(1), rate_kp(2))
%     fprintf('R1c %.8f, R2c %.8f, ', rate_c(1), rate_c(2))
%     fprintf('C1c %.8f, C2c %.8f, ', common_rates(1), common_rates(2));
%     fprintf('Rov = %.8f, ', u_k'*(rate_kp+common_rates));
%     fprintf('||Pc|| %.8f, ||Pk|| %.8f\n',trace(p_c*p_c'), trace(p_k*p_k'))
%     
    CVXstruct_SCA_ref(seed_value).iter = seed_value;
    CVXstruct_SCA_ref(seed_value).Rp = rate_kp;
    CVXstruct_SCA_ref(seed_value).Rc = rate_c;
    CVXstruct_SCA_ref(seed_value).Cc = common_rates;
    CVXstruct_SCA_ref(seed_value).Rov = u_k'*(rate_kp+common_rates);
    CVXstruct_SCA_ref(seed_value).P = [trace(p_c*p_c'), trace(p_k*p_k')];
    
    %% WMMSE
    
    [opt_val, p_k, p_c, common_rates] = ...
        CVX_opt_rsma_parameters_WMMSE(p_k_IC, p_c_IC,...
        K, N_T, u_k, Rth, Pt, h_ov_k, varianceNoise);
    
%     fprintf('WMMSE ref: ')
%     fprintf('opt value: %2.4f at %3d, ', opt_val, jj);
    [rate_c, rate_kp] = compute_rates(s, h_T_U_PL, h_R_U_PL, G, K, varianceNoise, p_c, p_k);
    
%     fprintf('R1p %.8f, R2p %.8f, ', rate_kp(1), rate_kp(2))
%     fprintf('R1c %.8f, R2c %.8f, ', rate_c(1), rate_c(2))
%     fprintf('C1c %.8f, C2c %.8f, ', common_rates(1), common_rates(2));
%     fprintf('Rov = %.8f, ', u_k'*(rate_kp+common_rates));
%     fprintf('||Pc|| %.8f, ||Pk|| %.8f\n',trace(p_c*p_c'), trace(p_k*p_k'))
    
    CVXstruct_WMMSE(seed_value).iter = seed_value;
    CVXstruct_WMMSE(seed_value).Rp = rate_kp;
    CVXstruct_WMMSE(seed_value).Rc = rate_c;
    CVXstruct_WMMSE(seed_value).Cc = common_rates;
    CVXstruct_WMMSE(seed_value).Rov = u_k'*(rate_kp+common_rates);
    CVXstruct_WMMSE(seed_value).P = [trace(p_c*p_c'), trace(p_k*p_k')];
end
name = 'rsma_cvx_comparison';
save(sprintf('./workspace_seed_list_%d_%s.mat', N, name))
%%
figure;
plot(1:N,[CVXstruct_SCA.Rov], '--x'); hold on;
plot(1:N,[CVXstruct_SCA_ref.Rov], '--o');
plot(1:N,[CVXstruct_WMMSE.Rov], '--d');
axis square; grid on;
leg = legend('SCA', 'SCA ref', 'WMMSE')

%%
