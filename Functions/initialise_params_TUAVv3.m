%%
% seed_value = 1;
rng(seed_value);
cpt = get(groot,'defaultLineLineWidth');
if cpt < 1
    set(groot,'defaultLineLineWidth',1.0)
end

% disp('########')
%%

% Based on: A Downlink Coverage Scheme of Tethered UAV (Table 2). 

% K = 2; % Number of users
N_T = 4; % Number of transmit antennas
N_UE = 1; % Number of antennas for each user;

Pt = 20*N_T; % 20 W per antenna.
f = 3.5e9; lambda = 3e8/f;

Thether_length = 100; % m

% N_R = 128; % Number of RE.

tau = 1; % amount of power split between commmon (1-tau)*Pt and private tau*Pt
% Weighted Sum-Rate Maximization for Rate-Splitting Multiple Access Based 
% Secure Communication

%% TUAV position
h_B = 30; % m
q_TUAV = [0, 0, h_B+Thether_length]';
% cable can be up to 100 m

%% RIS position

q_RIS = [10, 0, 30]';

elevationAngleRIS = atan((q_TUAV(3) - q_RIS(3))/norm(q_TUAV(1:2)-q_RIS(1:2)));

%% UEs position

x_pos = 5*(rand(K, 1)*2);
y_pos = 5*(rand(K, 1)*2-1);
z_pos = 1.5*ones(K, 1);
q_UEs = [x_pos, y_pos, z_pos]';

elevationAngleUEs = atan((q_TUAV(3) - q_UEs(3, :))./sqrt(sum((q_TUAV(1:2)-q_UEs(1:2, :)).^2)));

%%
deltaTxt = 0.3;
if pp
    figure;
    plot3(q_TUAV(1), q_TUAV(2), q_TUAV(3), '*'); hold on;
    plot3(q_RIS(1), q_RIS(2), q_RIS(3), 'ro')
    plot3(q_UEs(1, :), q_UEs(2, :), q_UEs(3, :),'ks')
    for k=1:K
        text(q_UEs(1, k)+deltaTxt,q_UEs(2, k)+deltaTxt,sprintf('%d', k))
    end
    grid on
    legend({'TUAV', 'RIS', 'UEs'})
end

%% distances
% direct path
d_T_U = sqrt(sum(q_UEs.^2 + q_TUAV.^2, 1));
% TUAV - RIS
d_T_R = sqrt(sum(q_RIS.^2 + q_TUAV.^2, 1));
% RIS - UE
d_R_U = sqrt(sum(q_UEs.^2 + q_RIS.^2, 1));


%% TUAV - UEs channel (direct)
urban = [9.61, 0.16];
denseUrban = [12.08, 0.11];
highriseUrban = [27.23, 0.08];
probabilityLOS = @(x, elevationAngle) 1./(1 + x(1)*exp(-x(2)*(elevationAngle - x(1))));

a1 = 2-3.5;
b1 = 3.5;

alphaTUAV2UE = a1*probabilityLOS(highriseUrban, elevationAngleUEs) + b1;

%%

% alphaTUAV2UE = 3.5;

KNLOS_factor_db = -10; % k factor depends on the position of the user
% a3 = 5; % 5dB
% b3 = 2/pi * (log(15/a3));
% KfactorTUAV2UE = a3 * exp(b3 * elevationAngleUEs)
KNLOS_factor = 10.^(KNLOS_factor_db/10);

h_T_U_LOS_factor = 1; % 
h_T_U_NLOS_factor = complex(randn(N_T, K),randn(N_T, K))/sqrt(2);

h_T_U_NLOS = (sqrt(KNLOS_factor/(1+KNLOS_factor))*h_T_U_LOS_factor +...
    sqrt(1/(1+KNLOS_factor))*h_T_U_NLOS_factor);

h_T_U_PL = h_T_U_NLOS.*sqrt((lambda/4/pi)^2 * d_T_U.^(-alphaTUAV2UE));

%% TUAV - RIS channel

alpha_T_R = 2;
KLOS_factor_db = 10;
KLOS_factor = 10.^(KLOS_factor_db/10);

h_T_R_LOS_factor = 1; 
h_T_R_NLOS_factor = complex(randn(N_R, N_T),randn(N_R, N_T))/sqrt(2);

h_T_R_LOS = (sqrt(KLOS_factor./(1+KLOS_factor)).*h_T_R_LOS_factor +...
    sqrt(1./(1+KLOS_factor)).*h_T_R_NLOS_factor);


G = h_T_R_LOS.*sqrt((lambda/4/pi)^2 * d_T_R.^(-alpha_T_R));

%% RIS - UEs channel

alpha_R_U = 2.8;

KLOS_factor_db = 10;
KLOS_factor = 10.^(KLOS_factor_db/10);

h_R_U_LOS_factor = 1; % 
h_R_U_NLOS_factor = complex(randn(N_R, N_UE),randn(N_R, N_UE))/sqrt(2);

h_R_U_LOS = (sqrt(KLOS_factor./(1+KLOS_factor)).*h_R_U_LOS_factor +...
                sqrt(1./(1+KLOS_factor)).*h_R_U_NLOS_factor); 

h_R_U_PL = h_R_U_LOS.*sqrt((lambda/4/pi)^2 * d_R_U.^(-alpha_R_U));

%% RIS
% Review Initial condition for phase shifts
phi = randi([-180, 180], N_R, 1)*pi/180;
% phi = pi*ones(N_R, 1);
s = exp(1i*phi);
Theta = diag(s);

%% Noise Variance

B = 20e6; % 20 MHz of BW; review the parameters for urban area
NF = 7;
ThermalNoise = db2pow(-173.9 + NF)/1000;  %W per Hz
% -173.9 is the noise power at 20 degrees
% - 204 is the noise power at 0 kelvin degrees
varianceNoise = ThermalNoise * B;

%% Overall channel
% amplification of the parameters for numerical resolution purpose.
h_T_U_PL = h_T_U_PL * 1e6;
h_R_U_PL = h_R_U_PL * 1e6;
varianceNoise = varianceNoise * (1e6)^2;

h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';

%% Precoder

p_c_IC = AMBF_common_precoder(h_ov_k, Pt, tau);
% p_c_IC = SVD_common_precoder(h_ov_k, Pt, tau);
p_k_IC = RZF_private_precoder_matrix(h_ov_k, Pt, K, tau, N_T);
% p_k_IC = MRT_precoder_private_matrix(h_ov_k, Pt, K, tau);

P = [p_c_IC, p_k_IC];
assert(round(trace(P*P')) == Pt, 'power out of bounds')



%%

[rate_c, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
common_rates = min(rate_c)/K*ones(size(rate_c));
% [~, UEk] = max(rate_kp);
% loop = 1;
% cpt_old = 1e4;
% maxiter = 1000;
% iter = 0;
% while loop
%     iter = iter + 1;
%     H_tilde_k = h_R_U_PL(:, UEk);
%     U1 =  diag(H_tilde_k') * G; 
%     phi = -angle(U1*(h_T_U_PL(:, UEk) + U1'*s));
%     s = exp(1i*phi);
%     cpt = sum(abs(h_T_U_PL(:, UEk) + U1'*s).^2);
%     if abs(cpt-cpt_old) < 1e-4
%         loop = 0;
%     end
%     if iter > maxiter
%         break;
%     end
%     cpt_old = cpt;
%     Theta = diag(s);
%     h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
%     [~, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
%     [~, UEk] = min(rate_kp);
% end
% Theta = diag(s);
% h_ov_k = (h_T_U_PL' + h_R_U_PL' * Theta * G)';
% p_c_IC = AMBF_common_precoder(h_ov_k, Pt, tau);
% % p_c_IC = SVD_common_precoder(h_ov_k, Pt, tau);
% p_k_IC = RZF_private_precoder_matrix(h_ov_k, Pt, K, tau, N_T);
% % p_k_IC = MRT_precoder_private_matrix(h_ov_k, Pt, K, tau);
% [rate_c, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
% common_rates = min(rate_c)/K*ones(size(rate_c));
% Rov = sum(rate_kp + common_rates);