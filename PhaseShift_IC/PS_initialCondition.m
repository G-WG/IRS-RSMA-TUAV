function [s, iter, loop] = PS_initialCondition(s, h_T_U_PL, h_R_U_PL, G, K, N_T, Pt, tau, varianceNoise, UEselection)

% PS_initialCondition(s, h_T_U_PL, h_R_U_PL, G, K, N_T, Pt, tau, varianceNoise)
%
% This functions provides an initial condition for the phase shift of the
% REs in a MISO environment. It seeks the best IRS configuration based on
% one user. The user could be randomly selected or be selected by the one
% who has the higher or lower channel capacity. It is based on an iterative
% approach.
% 
% Input parameters:
% - s: Phase shift configuration. Size of [N_R, 1].
% - h_T_U_PL : It is the channel from the TUAV to the UEs, each colum 
% represent the channel of a user. Size of [N_T, K].
% - h_R_U_PL : It is the channel from the RIS to the UEs, each colum 
% represent the channel of a user. Size of [N_R, K].
% - G: It is the channel from the TUAV to the RIS. Size of [N_R, N_T].
% - K : Number of users
% - N_T: Number of Transmit antennas at the B
% - Pt : Transmit power in watts
% - tau : Portion of the transmit power dedicated to the private stream.
% Note that (1 - tau) is for the common stream. Therefore should be between
% 0 and 1.
% - varianceNoise: Noise variance in Watts.
%
% Output parameters:
% - p_k: Precoder matrix for the private stream
%
% A typical entry could be:
% - RZF_private_precoder_matrix(h_ov_k, 5, 2, 0.5, 5)
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/07/10
%
% Reference:
% - 


%% Initial conditions
Theta = diag(s);
h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
p_c_IC = AMBF_common_precoder(h_ov_k, Pt, tau);
p_k_IC = RZF_private_precoder_matrix(h_ov_k, Pt, K, tau, N_T);
[~, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
% common_rates = min(rate_c)/K*ones(size(rate_c));

%% User selection
if ~isnumeric(UEselection)
    if all(UEselection == "max")
        [~, UEk] = max(rate_kp);
    elseif all(UEselection == 'min')
        [~, UEk] = min(rate_kp);
    elseif all(UEselection == 'random')
        UEk = randi([1, K]);
    end
else
    if UEselection <= K
        UEk = UEselection;
    else
        UEk = randi([1, K]);
    end
end

%% Iteration
loop = 1;
cpt_old = 1e4;
maxiter = 1000;
iter = 0;

while loop
    iter = iter + 1;
    H_tilde_k = h_R_U_PL(:, UEk);
    U1 =  diag(H_tilde_k') * G; 
    phi = -angle(U1*(h_T_U_PL(:, UEk) + U1'*s));
    s = exp(1i*phi);
    cpt = sum(abs(h_T_U_PL(:, UEk) + U1'*s).^2);
    if abs(cpt-cpt_old) < 1e-4
        loop = 0;
    end
    if iter > maxiter
        break;
    end
    cpt_old = cpt;
    Theta = diag(s);
    h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
    [~, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
    [~, UEk] = min(rate_kp);
end

end