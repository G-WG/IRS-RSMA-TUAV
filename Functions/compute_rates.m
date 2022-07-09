function [rate_c, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c, p_k)

% Theta = diag(s);
% h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
rate_kp = zeros(K, 1);
for k = 1:K
    interference_term = 0;
    for jj = 1:K
        if k ~= jj
            interference_term = interference_term + ...
            power(abs(h_ov_k(:, k)'*p_k(:, jj)), 2);
        end
    end
    rate_kp(k) = log2(1 + abs(h_ov_k(:, k)'*p_k(:, k))^2/(interference_term + varianceNoise));
end

rate_c = zeros(K, 1);
for k = 1:K
    rate_c(k) = log2(1 + (abs(h_ov_k(:, k)' * p_c).^2)/(sum(abs(h_ov_k(:, k)' * p_k).^2) + varianceNoise));
end

% end function
end