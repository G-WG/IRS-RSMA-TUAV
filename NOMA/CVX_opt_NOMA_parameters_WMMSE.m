function [opt_val, precoder_matrix_k, precoder_c, common_rates] = ...
        CVX_opt_NOMA_parameters_WMMSE(precoder_matrix_k_IC, precoder_c_IC,...
        K, N_T, u_k, Rth, Pt, h_ov_k, varianceNoise, UE_NOMA_common, UE_NOMA_private)
%% Dubts
% Needs to reformulate CVX opt, due to the interference. One user receive
% s12 and decode it treating the rest as noise, and the other uses SIC to
% remove that message and decode s2 (for instance) without interference


%% Initial lambda and e MMSE optimal conditions

lambda_c_IC = zeros(K, 1);
lambda_p_IC = zeros(K, 1);
e_c_IC = zeros(K, 1);
e_p_IC = zeros(K, 1);
T_k_c = zeros(K, 1);
I_k_c = zeros(K, 1);
T_k_p = zeros(K, 1);
I_k_p = zeros(K, 1);

for k = 1:K
    interference_term = 0;
    for jj = 1:K
        if k ~= jj && jj ~= UE_NOMA_common
            interference_term = interference_term + ...
            power(abs(h_ov_k(:, k)'*precoder_matrix_k_IC(:, jj)), 2);
        end
    end

    % averaged received power at each user-k, equation (40)
    I_k_p(k) = interference_term + varianceNoise;
    if k == UE_NOMA_private
        I_k_c(k) = I_k_p(k) + pow_abs(h_ov_k(:, k)'*precoder_matrix_k_IC(:, k), 2);
    else
        I_k_c(k) = I_k_p(k);
    end
    T_k_p(k) = I_k_c(k);
    T_k_c(k) = T_k_p(k) + pow_abs(h_ov_k(:, k)'*precoder_c_IC, 2);

    % MMSE equalisers, equations (43) and (44)
    e_c_IC(k) = precoder_c_IC'*h_ov_k(:, k) / T_k_c(k);
    e_p_IC(k) = precoder_matrix_k_IC(:, k)'*h_ov_k(:, k) / T_k_p(k);

    % Weights of the augmented WMSE of the common and private stream,
    % equation (51) with value of eq. (45) and (46)
    lambda_c_IC(k) = T_k_c(k)/I_k_c(k);
    lambda_p_IC(k) = T_k_p(k)/I_k_p(k);
    
end


%% CVX

cvx_begin quiet

    % variables
    variable v_common_rates(K)
    variable precoder_matrix_k(N_T, K) complex
    variable precoder_c(N_T, 1) complex

    T_k_c = cvx(zeros(K, 1));
    I_k_c = cvx(zeros(K, 1));
    T_k_p = cvx(zeros(K, 1));
    I_k_p = cvx(zeros(K, 1));

    varepsilon_k_c = cvx(zeros(K, 1));
    varepsilon_k_p = cvx(zeros(K, 1));

    xi_k_c = cvx(zeros(K, 1));
    xi_k_p = cvx(zeros(K, 1));

    for k = 1:K
        
        interference_term = 0;
        for jj = 1:K
            if k ~= jj && jj ~= UE_NOMA_common
                interference_term = interference_term + ...
                square_abs(h_ov_k(:, k)'*precoder_matrix_k(:, jj));
            end
        end
        I_k_p(k) = interference_term + varianceNoise;
        if k == UE_NOMA_private
            I_k_c(k) = I_k_p(k) + square_abs(h_ov_k(:, k)'*precoder_matrix_k(:, k));
        else
            I_k_c(k) = I_k_p(k);
        end
        T_k_p(k) = I_k_c(k);
        T_k_c(k) = T_k_p(k) + square_abs(h_ov_k(:, k)'*precoder_c);

        % MSE calculation, equationa (41) and (42)
        varepsilon_k_c(k) = square_abs(e_c_IC(k))*T_k_c(k) - ...
                            2*real(e_c_IC(k)*h_ov_k(:, k)'*precoder_c) + 1;
        varepsilon_k_p(k) = square_abs(e_p_IC(k))*T_k_p(k) - ...
                            2*real(e_p_IC(k)*h_ov_k(:, k)'*precoder_matrix_k(:, k)) + 1;

        % Augmented WMSE equation (49)
        xi_k_c(k) = lambda_c_IC(k)*varepsilon_k_c(k) - log2(lambda_c_IC(k));
        if k == UE_NOMA_private
            xi_k_p(k) = lambda_p_IC(k)*varepsilon_k_p(k) - log2(lambda_p_IC(k));
        end
        
    end

    objectFunction = u_k'*[v_common_rates + xi_k_p];

    minimize( objectFunction )

    subject to
        % constraint (a)
        sum(v_common_rates) + 1 >= xi_k_c;

        % constraint (b)
        xi_k_p + v_common_rates <= 1 - Rth;
%         if UE_NOMA_private == 2
%             [v_common_rates(UE_NOMA_common); xi_k_p(UE_NOMA_private)] <= 1 - Rth;
%         else
%             [xi_k_p(UE_NOMA_private); v_common_rates(UE_NOMA_common)] <= 1 - Rth;
%         end

        % constraint (c)

        P = [precoder_c, precoder_matrix_k]; 
        sum(sum_square_abs(P)) <= Pt; % trace(p * p') <= Pt;
        v_common_rates <= 0;

cvx_end

opt_val = cvx_optval;
common_rates = -v_common_rates;

end