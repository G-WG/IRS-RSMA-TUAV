function p_k = MRT_precoder_private_matrix(h_ov_k, Pt, K, tau)

% MRT_precoder_private_matrix(h_ov_k, Pt, K, tau)
%
% This functions provides the precoder matrix for the private streams when
% RSMA scheme is used. It is based on MRT.
% 
% Input parameters:
% - h_ov_k : It is the overall channel, each colum represent the channel of
% a user.
% - Pt : Transmit power in watts
% - K : Number of users
% - tau : Portion of the transmit power dedicated to the private stream.
% Note that (1 - tau) is for the common stream. Therefore should be between
% 0 and 1.
%
% Output parameters:
% - p_k: Precoder matrix for the private stream
%
% A typical entry could be:
% - MRT_precoder_private_matrix(h_ov_k, 5, 2, 0.5)
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/06/21
%
% Reference:
% - Mao, Y., Clerckx, B. & Li, V.O. Rate-splitting multiple access for 
%   downlink communication systems: bridging, generalizing, and 
%   outperforming SDMA and NOMA. J Wireless Com Network 2018, 133 (2018). 
%   https://doi.org/10.1186/s13638-018-1104-7

p_k = h_ov_k./sqrt(sum(abs(h_ov_k.^2))) * sqrt(tau*Pt/(K));
assert(round(trace(p_k*p_k'))==round(tau*Pt), 'Power budget error')

end