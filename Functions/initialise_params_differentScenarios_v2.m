%%
% seed_value = 5; pp = 1; K = 2; N_R = 128; tau = 0.8; precoderIC = 1;
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

% tau = 0.8; % amount of power split between commmon (1-tau)*Pt and private tau*Pt
% Weighted Sum-Rate Maximization for Rate-Splitting Multiple Access Based 
% Secure Communication

wB = 20; % distance between building centers
wS = wB/2; % street width
widthBuilding = wB-wS;
cornerBuildings = [widthBuilding/2, widthBuilding/2; -widthBuilding/2, widthBuilding/2];
%% RIS position

q_RIS = [widthBuilding/2 + wS, 0, 20]';

elevationAngleRIS = atan((q_TUAV(3) - q_RIS(3))/norm(q_TUAV(1:2)-q_RIS(1:2)));

%% UEs position

q_UEs = zeros(3, K);
q_UEs(:, 1) = [q_RIS(1)-3, q_RIS(2)+3, 1.5];
q_UEs(:, 2) = [q_RIS(1)-wS+1, q_RIS(2)-1, 1.5];
% q_UEs(:, 3) = [q_RIS(1)-wS-widthBuilding-2, q_RIS(2)-3, 1.5];
elevationAngleUEs = atan((q_TUAV(3) - q_UEs(3, :))./sqrt(sum((q_TUAV(1:2)-q_UEs(1:2, :)).^2)));

q = zeros(2, K); q(1, :) = widthBuilding/2; q(2, :) = q_UEs(2, :);
elevationAngleUEs2BuildingCenter = atan((h_B - q_UEs(3, :))./sqrt(sum(q-q_UEs(1:2, :)).^2));

% elevationAngleUEs2BuildingCenter < elevationAngleUEs
q = zeros(2, 1); q(1, :) = widthBuilding/2; q(2, :) = q_RIS(2, :);
elevationAngleRIS2BuildingCenter = atan((h_B - q_RIS(3, :))./sqrt(sum(q-q_RIS(1:2, :)).^2));

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
    legend({'TUAV', 'RIS', 'UEs'})

    r = rectangle('Position',[-widthBuilding/2 -widthBuilding/2 widthBuilding widthBuilding]);
    rectangle('Position',[-widthBuilding/2+wB -widthBuilding/2 widthBuilding widthBuilding]);
    xlim([-15, 30])
    ylim([-15, 15])
end

%% Define visibility cones

thetaRightCorner2RIS = atan((q_RIS(2) - cornerBuildings(1, 1))/(q_RIS(1) - cornerBuildings(1, 2)));
thetaLeftCorner2RIS = atan((q_RIS(2) - cornerBuildings(2, 1))/(q_RIS(1) - cornerBuildings(2, 2)));
thetaTUAV2RIS = atan((q_RIS(2) - q_TUAV(2))/(q_RIS(1) - q_TUAV(1)));

if thetaTUAV2RIS < thetaRightCorner2RIS & thetaTUAV2RIS > thetaLeftCorner2RIS
    disp('TUAV is in the cone')
    % In order to know whether the TUAV and RIS has LoS, elevation angles
    % must be compared
else
    disp('LOS TUAV-RIS')
end




