function [opt_val, precoder_matrix_k] = ...
        CVX_opt_SDMA_parameters_WMMSE(precoder_matrix_k_IC,...
        K, N_T, u_k, Rth, Pt, h_ov_k, varianceNoise)



%% Initial lambda and e MMSE optimal conditions


lambda_p_IC = zeros(K, 1);
e_p_IC = zeros(K, 1);
T_k_p = zeros(K, 1);
I_k_p = zeros(K, 1);

for k = 1:K
    interference_term = 0;
    for jj = 1:K
        if k ~= jj
            interference_term = interference_term + ...
            power(abs(h_ov_k(:, k)'*precoder_matrix_k_IC(:, jj)), 2);
        end
    end

    % averaged received power at each user-k, equation (40)
    I_k_p(k) = interference_term + varianceNoise;
    T_k_p(k) = I_k_p(k) + pow_abs(h_ov_k(:, k)'*precoder_matrix_k_IC(:, k), 2);

    % MMSE equalisers, equations (43) and (44)
    e_p_IC(k) = precoder_matrix_k_IC(:, k)'*h_ov_k(:, k) / T_k_p(k);

    % Weights of the augmented WMSE of the common and private stream,
    % equation (51) with value of eq. (45) and (46)
    lambda_p_IC(k) = T_k_p(k)/I_k_p(k);
    
end


%% CVX

cvx_begin quiet

    % variables
    variable precoder_matrix_k(N_T, K) complex

    T_k_p = cvx(zeros(K, 1));
    I_k_p = cvx(zeros(K, 1));

    varepsilon_k_p = cvx(zeros(K, 1));

    xi_k_p = cvx(zeros(K, 1));

    for k = 1:K
        
        interference_term = 0;
        for jj = 1:K
            if k ~= jj
                interference_term = interference_term + ...
                square_abs(h_ov_k(:, k)'*precoder_matrix_k(:, jj));
            end
        end
        I_k_p(k) = interference_term + varianceNoise;
        T_k_p(k) = I_k_p(k) + square_abs(h_ov_k(:, k)'*precoder_matrix_k(:, k));
        

        % MSE calculation, equationa (41) and (42)
        varepsilon_k_p(k) = square_abs(e_p_IC(k))*T_k_p(k) - ...
                            2*real(e_p_IC(k)*h_ov_k(:, k)'*precoder_matrix_k(:, k)) + 1;

        % Augmented WMSE equation (49)
        xi_k_p(k) = lambda_p_IC(k)*varepsilon_k_p(k) - log2(lambda_p_IC(k));
        
    end

    objectFunction = u_k'*(xi_k_p);

    minimize( objectFunction )

    subject to
        % constraint (a)
       
        % constraint (b)
        xi_k_p <= 1 - Rth;

        % constraint (c)

        P = [precoder_matrix_k]; 
        sum(sum_square_abs(P)) <= Pt; % trace(p * p') <= Pt;

cvx_end

opt_val = cvx_optval;

end