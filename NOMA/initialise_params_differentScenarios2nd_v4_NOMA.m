%%
% seed_value = 150; pp = 1; K = 2; N_R = 128; tau = 0.5; precoderIC = 1; h_B = 30; q_TUAV = [-5, 0, h_B+20]';
rng(seed_value);

%% INITIAL PARAMS

UE_NOMA_common = 2;
UE_NOMA_private = 1;
assert(UE_NOMA_common ~= UE_NOMA_private, 'Select different users');


% K = 2; % Number of users
N_T = 4; % Number of transmit antennas
N_UE = 1; % Number of antennas for each user;

% Pt = db2pow(Pt_dB)/1000;
Pt = 20*N_T; % 20 W per antenna.
f = 3.5e9; lambda = 3e8/f;

Thether_length = 100; % m

% N_R = 128; % Number of RE.

% tau = 0.8; % amount of power split between commmon (1-tau)*Pt and private tau*Pt
% Weighted Sum-Rate Maximization for Rate-Splitting Multiple Access Based 
% Secure Communication

stretchingFactor = 5;
wB = 20; % distance between building centers
wS = wB/2; % street width
widthBuilding = wB-wS;
% cornerBuildings = [widthBuilding/2+stretchingFactor, widthBuilding/2; -widthBuilding/2-stretchingFactor, -widthBuilding/2];
% cornerBuilding = [y, x; y, x]

xiList = 0:2:4;
yiList = -4:2:4;

cpt = 1;
for xi = xiList
    for yi = yiList
        cBx = -widthBuilding/2-xi*widthBuilding;
        cBy = -widthBuilding/2-stretchingFactor-yi*(widthBuilding+stretchingFactor);
        S(cpt).cornerBuildings = [cBy+widthBuilding+stretchingFactor*2, cBx+widthBuilding; ...
                                cBy, cBx];
        S(cpt).height = h_B;
        cpt = cpt + 1;
    end
end
%% RIS position

q_RIS = [widthBuilding/2 + wS, 0, 20]';

RIS2TUAVElevation = atan((q_TUAV(3) - q_RIS(3))/norm(q_TUAV(1:2)-q_RIS(1:2)));

%% UEs position

q_UEs = zeros(3, K);
q_UEs(:, 1) = [wS/2+1, q_RIS(2)+1, 1.5];
q_UEs(:, 2) = [-wS/2-1, q_RIS(2)-1, 1.5];
% q_UEs(:, 3) = [q_RIS(1)-wS-widthBuilding-2, q_RIS(2)-3, 1.5];
UEs2TUAVElevation = atan2((q_TUAV(3) - q_UEs(3, :)), sqrt(sum((q_TUAV(1:2)-q_UEs(1:2, :)).^2)));

q = zeros(2, K); q(1, :) = widthBuilding/2; q(2, :) = q_UEs(2, :);
elevationAngleUEs2BuildingCenter = atan((h_B - q_UEs(3, :))./sqrt(sum(q-q_UEs(1:2, :)).^2));


%%
deltaTxt = 0.5;

if pp
    figure;
    plot3(q_TUAV(1), q_TUAV(2), q_TUAV(3), '*'); hold on;
    plot3(q_RIS(1), q_RIS(2), q_RIS(3), 'rs')
    plot3(q_UEs(1, :), q_UEs(2, :), q_UEs(3, :),'ko')
    for k=1:K
        text(q_UEs(1, k)+deltaTxt,q_UEs(2, k)+deltaTxt,sprintf('%d', k))
    end
    grid on
    leg = legend({'TUAV', 'RIS', 'UEs'});

%     rectangle('Position',[-widthBuilding/2 -widthBuilding/2-stretchingFactor widthBuilding widthBuilding+stretchingFactor*2]);
%     rectangle('Position',[-widthBuilding/2+wB -widthBuilding/2-stretchingFactor widthBuilding widthBuilding+stretchingFactor*2]);
    cpt = 1;
    for xi = [-2, xiList]
        for yi = yiList
            rectangle('Position',[-widthBuilding/2-xi*widthBuilding ...
                -widthBuilding/2-stretchingFactor-yi*(widthBuilding+stretchingFactor) ...
                widthBuilding widthBuilding+stretchingFactor*2]);
            xc = -widthBuilding/2-xi*widthBuilding;
            yc = -widthBuilding/2-stretchingFactor-yi*(widthBuilding+stretchingFactor);
            if xi > -2
                text(xc+deltaTxt,yc+deltaTxt,sprintf('%d', cpt))
                cpt = cpt+1;
            end
        end
    end
    xlim([-100, 100])
    ylim([-100, 100])
    xlabel('m')
    ylabel('m')
    legbool=1;setPlotParams;
    axis square
end


%% RIS LoS check
objectPosition = q_UEs(:, 2) ;
object2TUAVElevation = UEs2TUAVElevation(2);
idxLOSBoleanArrayUE2 = isLOS(S, objectPosition, object2TUAVElevation, q_TUAV);
isUE2LOS = all(boolean(idxLOSBoleanArrayUE2));

objectPosition = q_UEs(:, 1) ;
object2TUAVElevation = UEs2TUAVElevation(1);
idxLOSBoleanArrayUE1 = isLOS(S, objectPosition, object2TUAVElevation, q_TUAV);
isUE1LOS = all(boolean(idxLOSBoleanArrayUE1));

objectPosition = q_RIS ;
object2TUAVElevation = RIS2TUAVElevation;
idxLOSBoleanArrayRIS = isLOS(S, objectPosition, object2TUAVElevation, q_TUAV);
isRISLOS = all(boolean(idxLOSBoleanArrayRIS));

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

alphaTUAV2UE = a1*probabilityLOS(highriseUrban, UEs2TUAVElevation) + b1;

% alphaTUAV2UE = 3.5;
a3 = 5; % 5dB
b3 = 2/pi * (log(15/a3));
KNLOS_factor_db = a3 * exp(b3 * UEs2TUAVElevation);

if ~isUE1LOS
    KNLOS_factor_db(1) = -100;
end
if ~isUE2LOS
    KNLOS_factor_db(2) = -100;
end

KNLOS_factor = 10.^(KNLOS_factor_db/10);

h_T_U_LOS_factor = 1; % 
h_T_U_NLOS_factor = complex(randn(N_T, K),randn(N_T, K))/sqrt(2)/sqrt(N_T*K);

h_T_U_NLOS = (sqrt(KNLOS_factor./(1+KNLOS_factor)).*h_T_U_LOS_factor +...
    sqrt(1./(1+KNLOS_factor)).*h_T_U_NLOS_factor);

h_T_U_PL = h_T_U_NLOS.*sqrt((lambda/4/pi)^2 * d_T_U.^(-alphaTUAV2UE));

%% TUAV - RIS channel

alpha_T_R = 2;


KLOS_factor_db = 10;
if ~isRISLOS
    KLOS_factor_db = -100;
end
KLOS_factor = 10.^(KLOS_factor_db/10);

h_T_R_LOS_factor = 1; 
h_T_R_NLOS_factor = complex(randn(N_R, N_T),randn(N_R, N_T))/sqrt(2)/sqrt(N_R*N_T);

h_T_R_LOS = (sqrt(KLOS_factor./(1+KLOS_factor)).*h_T_R_LOS_factor +...
    sqrt(1./(1+KLOS_factor)).*h_T_R_NLOS_factor);


G = h_T_R_LOS.*sqrt((lambda/4/pi)^2 * d_T_R.^(-alpha_T_R));

%% RIS - UEs channel

alpha_R_U = 2.8;

KLOS_factor_db = [10 -100];
KLOS_factor = 10.^(KLOS_factor_db/10);

h_R_U_LOS_factor = 1; % 
h_R_U_NLOS_factor = complex(randn(N_R, K),randn(N_R, K))/sqrt(2)/sqrt(N_R*K);

h_R_U_LOS = (sqrt(KLOS_factor./(1+KLOS_factor)).*h_R_U_LOS_factor +...
                sqrt(1./(1+KLOS_factor)).*h_R_U_NLOS_factor); 

h_R_U_PL = h_R_U_LOS.*sqrt((lambda/4/pi)^2 * d_R_U.^(-alpha_R_U));

%%

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
if precoderIC == 1
    p_c_IC = AMBF_common_precoder(h_ov_k, Pt, tau);
    p_k_IC = RZF_private_precoder_matrix(h_ov_k, Pt, K, tau, N_T);
    p_k_IC(:, UE_NOMA_common) = 0;
    p_k_IC = p_k_IC/norm(p_k_IC) * sqrt(tau*Pt);
    
else
    p_c_IC = SVD_common_precoder(h_ov_k, Pt, tau);
    p_k_IC = MRT_precoder_private_matrix(h_ov_k, Pt, K, tau);
    p_k_IC(:, UE_NOMA_common) = 0;
    p_k_IC = p_k_IC/norm(p_k_IC) * sqrt(tau*Pt);
end
P = [p_c_IC, p_k_IC];
assert(round(trace(P*P')) == Pt, 'power out of bounds')



%%

[rate_c, rate_kp] = compute_rates(h_ov_k, K, varianceNoise, p_c_IC, p_k_IC);
% common_rates = zeros(size(rate_c)); common_rates(UE_NOMA_common) = rate_c(UE_NOMA_common);
common_rates = min(rate_c)*ones(K, 1)/K;

wsr = sum(rate_kp + common_rates);

