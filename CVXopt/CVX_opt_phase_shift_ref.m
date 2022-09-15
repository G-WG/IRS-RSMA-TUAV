function [opt_val, s, eta_p, beta_p] = CVX_opt_phase_shift_ref(s_IC, ...
    eta_p_IC, beta_p_IC, C, K, N_R, h_T_U_PL, h_R_U_PL, G, p_c, p_k, ...
    varianceNoise, Lambda_c, Lambda_p)


% [opt_val, s, eta_p, beta_p, beta_c] = CVX_opt_phase_shift(s_IC, ...
%    eta_IC, beta_IC, C, K, h_T_U_PL, h_R_U_PL, G, p_c, p_k, varianceNoise)
%
% This functions performs a phase shift optimsation using CVX.
% 
% Input parameters:
% - s_IC : Initial condition for the phase shift. Size of [N_R, 1]
% - eta_p_IC : Initial condition for the slack variable eta. Size of [K, 1]
% - beta_p_IC : Initial condition for the slack variable beta. Size of [K, 1]
% - C : Penalty constant
% - K :  number of users.
% - h_T_U_PL : Channel from the TUAV to the user.
% - h_R_U_PL : Channel from the IRS to the user.
% - G : Channel from the TUAV to the IRS.
% - p_c : Precoder vector for the common stream.
% - p_k : Precoder matrix for the private stream. Size of [N_T, K], N_T
% number of Transmit antennas.
% - varianceNoise : Noise variance in Watts.
%
% Output parameters:
% - opt_val: Optimal value from the CVX optimisation problem
% - s : optimal phase shift.
% - eta_p : optimal eta_p.
% - beta_c : optimal beta_c.
% - beta_p : optimal beta_p.
%
% A typical entry could be:
% - [opt_val, s, eta_p, beta_p, beta_c] = CVX_opt_phase_shift(s_IC, ...
%    eta_IC, beta_IC, C, K, h_T_U_PL, h_R_U_PL, G, p_c, p_k, varianceNoise)
% with:
%   - s_IC = exp(1i*randi([-180, 180], N_R, 1))
%   - eta_IC = zeros(N_R)
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/06/21


%%
s_old = s_IC;
eta_p_old  = eta_p_IC;
beta_p_old = beta_p_IC;

cvx_begin quiet
    % cvx_precision high
    % variables
    variable s(N_R) complex
    variable eta_p(K) nonnegative
    variable beta_p(K) nonnegative

    % auxiliary variables
    w_obj_rate = sum(eta_p);
    w_obj_penalty = 2*real(s_old'*(s - s_old));
%     w_obj_penalty = sum(abs(s_old.^2) - 1 + real(2*conj(s_old).*(s - s_old)));

    % maximize
    maximize( w_obj_rate + C*w_obj_penalty)

    subject to
        % constraint (P4.a)
        abs(s) <= 1;

        for k = 1:K
            
            interference_term = 0;
            % common streams
            hdkc = h_T_U_PL(:, k)'*p_c;
            tkc = conj(diag(h_R_U_PL(:, k)')*G*p_c); % tkcs = tkc'*s;

            % private streams
            tk = diag(h_R_U_PL(:, k)')*G;
            tkk = conj(diag(h_R_U_PL(:, k)')*G*p_k(:, k)); % tkks = tkk'*s;
            hdkk = h_T_U_PL(:, k)'*p_k(:, k);

            for jj = 1:K
                if k ~= jj
                    interference_term = interference_term + ...
                        pow_abs( (h_T_U_PL(:, k)'*p_k(:, jj) +  ...
                        transpose(tk*p_k(:, jj))*s), 2 );
                end
            end
                % constraint (19)
            beta_p(k) >= interference_term + varianceNoise;
            % constraint (20)
            pow_abs( (hdkk + tkk'*s_old), 2 ) + ...
                2*real((hdkk + tkk'*s_old)'*tkk'*(s - s_old)) >= ...
                1/4*( power(beta_p(k) + eta_p(k), 2) - ...
                2*(beta_p_old(k) - eta_p_old(k))*(beta_p(k) - eta_p(k)) + ...
                power(beta_p_old(k) - eta_p_old(k), 2) );
                
            % constraint (21)
%             beta_c(k) >= beta_p(k) + pow_abs((hdkk + tkk'*s), 2);
            % contraint (23)
            pow_abs( (hdkc + tkc'*s_old), 2) + ...
                2*real((hdkc + tkc'*s_old)'*tkc'*(s - s_old)) >= ...
                    Lambda_c * (interference_term + varianceNoise + ...
                    pow_abs( (hdkk + tkk'*s), 2 ));
            
            
            
        end
        % constraint (P2.c)
        eta_p >= Lambda_p;
cvx_end

opt_val = cvx_optval;
% if opt_val == Inf || isnan(opt_val)
%     disp('inf-nan')
% end
