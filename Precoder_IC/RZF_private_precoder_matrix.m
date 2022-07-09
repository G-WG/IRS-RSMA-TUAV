function p_k = RZF_private_precoder_matrix(h_ov_k, Pt, K, tau, N_T)

% RZF_private_precoder_matrix(h_ov_k, Pt, K, tau, N_T)
%
% This functions provides the precoder matrix for the private streams when
% RSMA scheme is used. It is based on regularized Zero Forcing assuming
% perfect CSIT.
% 
% Input parameters:
% - h_ov_k : It is the overall channel, each colum represent the channel of
% a user.
% - Pt : Transmit power in watts
% - K : Number of users
% - tau : Portion of the transmit power dedicated to the private stream.
% Note that (1 - tau) is for the common stream. Therefore should be between
% 0 and 1.
% - N_T: Number of Transmit antennas at the BS.
%
% Output parameters:
% - p_k: Precoder matrix for the private stream
%
% A typical entry could be:
% - RZF_private_precoder_matrix(h_ov_k, 5, 2, 0.5, 5)
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/06/21
%
% Reference:
% - G. Lu, L. Li, H. Tian and F. Qian, "MMSE-Based Precoding for Rate 
%   Splitting Systems With Finite Feedback," in IEEE Communications Letters, 
%   vol. 22, no. 3, pp. 642-645, March 2018, doi: 10.1109/LCOMM.2017.2785221.

h_ov_k_norm = h_ov_k./sqrt(sum(abs(h_ov_k.^2)));
eta = sqrt(K/(norm((h_ov_k_norm*h_ov_k_norm' + K/((tau*Pt/K)*N_T)*eye(N_T))\h_ov_k_norm, "fro")^2));
p_k = (h_ov_k_norm*h_ov_k_norm' + K/((tau*Pt/K)*N_T)*eye(N_T))\h_ov_k_norm * eta * sqrt(tau*Pt/(K));
assert(round(trace(p_k*p_k'))==round(tau*Pt), 'Power budget error')

end