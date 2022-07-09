function p_c = SVD_common_precoder(h_ov_k, Pt, tau)

% SVD_common_precoder(h_ov_k, Pt, tau)
%
% This functions provides the precoder vector for the common stream when
% RSMA scheme is used. It is based on SVD decomposition.
% 
% Input parameters:
% - h_ov_k : It is the overall channel, each colum represent the channel of
% a user.
% - Pt : Transmit power in watts
% - tau : Portion of the transmit power dedicated to the private stream.
% Note that (1 - tau) is for the common stream. Therefore should be between
% 0 and 1.
%
% Output parameters:
% - p_c: Precoder vector for the common stream
%
% A typical entry could be:
% - SVD_common_precoder(h_ov_k, 5, 0.9)
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/06/21
%
% Reference:
% - Mao, Y., Clerckx, B. & Li, V.O. Rate-splitting multiple access for 
%   downlink communication systems: bridging, generalizing, and 
%   outperforming SDMA and NOMA. J Wireless Com Network 2018, 133 (2018). 
%   https://doi.org/10.1186/s13638-018-1104-7

[U,~,~] = svd(h_ov_k);
p_c = U(:, 1) * sqrt((1-tau)*Pt);
assert(round(trace(p_c*p_c'))==round((1-tau)*Pt), 'Power budget error')

end
