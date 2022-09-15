function [opt_val, precoder_c, common_rates] = ...
        CVX_opt_NOMA_parameters_WMMSE(precoder_c_IC,...
        K, N_T, u_k, Rth, Pt, h_ov_k, varianceNoise)

% NOMA:  Non-Orthogonal Multiple Access
% It is the common precoder of an RSMA scheme using 1L.

%% Initial lambda and e MMSE optimal conditions

lambda_c_IC = zeros(K, 1);
e_c_IC = zeros(K, 1);
T_k_c = zeros(K, 1);
I_k_c = zeros(K, 1);

for k = 1:K

    I_k_c(k) = varianceNoise;
    T_k_c(k) = pow_abs(h_ov_k(:, k)'*precoder_c_IC, 2);

    % MMSE equalisers, equations (43) and (44)
    e_c_IC(k) = precoder_c_IC'*h_ov_k(:, k) / T_k_c(k);

    % Weights of the augmented WMSE of the common and private stream,
    % equation (51) with value of eq. (45) and (46)
    lambda_c_IC(k) = T_k_c(k)/I_k_c(k);
    
end


%% CVX

cvx_begin quiet
%     cvx_precision high
    % variables
    variable v_common_rates(K)
    variable precoder_c(N_T, 1) complex

    T_k_c = cvx(zeros(K, 1));
    I_k_c = cvx(zeros(K, 1));

    varepsilon_k_c = cvx(zeros(K, 1));

    xi_k_c = cvx(zeros(K, 1));

    for k = 1:K
        
        I_k_c(k) = varianceNoise;
        T_k_c(k) = square_abs(h_ov_k(:, k)'*precoder_c);

        % MSE calculation, equationa (41) and (42)
        varepsilon_k_c(k) = square_abs(e_c_IC(k))*T_k_c(k) - ...
                            2*real(e_c_IC(k)*h_ov_k(:, k)'*precoder_c) + 1;

        % Augmented WMSE equation (49)
        xi_k_c(k) = lambda_c_IC(k)*varepsilon_k_c(k) - log2(lambda_c_IC(k));
        
    end

    objectFunction = u_k'*(v_common_rates);

    minimize( objectFunction )

    subject to
        % constraint (a)
        sum(v_common_rates) + 1 >= max(xi_k_c);

        % constraint (b)
        v_common_rates <= -Rth;

        % constraint (c)

        P = precoder_c; 
        sum(sum_square_abs(P)) <= Pt; % trace(p * p') <= Pt;
        v_common_rates <= 0;

cvx_end

opt_val = cvx_optval;
common_rates = -v_common_rates;

end