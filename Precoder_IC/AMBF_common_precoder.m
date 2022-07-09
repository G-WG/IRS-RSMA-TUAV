function p_c = AMBF_common_precoder(h_ov_k, Pt, tau)

% AMBF_common_precoder(h_ov_k, Pt, tau)
%
% This functions provides the precoder vector for the common stream when
% RSMA scheme is used. It is based on computing the Average Matched
% beamforming
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
% - MRT_precoder_private_matrix(h_ov_k, 5, 2, 0.9)
%
% Author: Maximiliano Rivera --  marivera3@uc.cl
% Version: v1.0 2022/06/21
%
% Reference:
% - G. Lu, L. Li, H. Tian and F. Qian, "MMSE-Based Precoding for Rate 
%   Splitting Systems With Finite Feedback," in IEEE Communications Letters, 
%   vol. 22, no. 3, pp. 642-645, March 2018, doi: 10.1109/LCOMM.2017.2785221.

p_c = sqrt((1-tau)*Pt) * sum(h_ov_k, 2)/norm(sum(h_ov_k, 2));
assert(round(trace(p_c*p_c'))==round((1-tau)*Pt), 'Power budget error')

end