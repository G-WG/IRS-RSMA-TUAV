%%
% seed_value = 5; pp = 1; K = 2; N_R = 128; tau = 0.8; precoderIC = 1;
rng(seed_value);
% cpt = get(groot,'defaultLineLineWidth');
% if cpt < 1
%     set(groot,'defaultLineLineWidth',1.0)
% end

% disp('########')
%%

% Based on: A Downlink Coverage Scheme of Tethered UAV (Table 2). 

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
cornerBuildings = [widthBuilding/2+stretchingFactor, widthBuilding/2; -widthBuilding/2-stretchingFactor, -widthBuilding/2];
% cornerBuilding = [y, x; y, x]

xiList = [0:2:4];
yiList = [-4:2:4];

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
    legend({'TUAV', 'RIS', 'UEs'})

%     rectangle('Position',[-widthBuilding/2 -widthBuilding/2-stretchingFactor widthBuilding widthBuilding+stretchingFactor*2]);
%     rectangle('Position',[-widthBuilding/2+wB -widthBuilding/2-stretchingFactor widthBuilding widthBuilding+stretchingFactor*2]);

    for xi = [-2, xiList]
        for yi = yiList
            rectangle('Position',[-widthBuilding/2-xi*widthBuilding ...
                -widthBuilding/2-stretchingFactor-yi*(widthBuilding+stretchingFactor) ...
                widthBuilding widthBuilding+stretchingFactor*2]);
        end
    end
    xlim([-100, 100])
    ylim([-100, 100])
    axis square
end

%%

idxRIS = ones(length(S), 1);

for ii = 1:length(S)
    if q_TUAV(2) < -wB 
        % The cone is built by the left upper corner and right lower corner of
        % the buildings
        if S(ii).cornerBuildings(1, 1) <= -wB
            thetaRightCorner2RIS = atan((q_RIS(2) - S(ii).cornerBuildings(1, 1))/(q_RIS(1) - S(ii).cornerBuildings(2, 2))); % y / x
            thetaLeftCorner2RIS = atan((q_RIS(2) - S(ii).cornerBuildings(2, 1))/(q_RIS(1) - S(ii).cornerBuildings(1, 2)));
            thetaTUAV2RIS = atan2((q_RIS(2) - q_TUAV(2)), (q_RIS(1) - q_TUAV(1)));
            if thetaTUAV2RIS > thetaRightCorner2RIS && thetaTUAV2RIS < thetaLeftCorner2RIS
            %     disp('TUAV is in the cone')
                % In order to know whether the TUAV and RIS has LoS, elevation angles
                % must be compared
                RIS2Building = (q_RIS(1) - S(ii).cornerBuildings(1, 2)) / cos(thetaTUAV2RIS);
                RIS2BuildingElevation = atan((S(ii).height-q_RIS(3))/RIS2Building);
                idxRIS(ii) = RIS2TUAVElevation >= RIS2BuildingElevation;
            end
        end
    
    elseif q_TUAV(2) > wB
        if S(ii).cornerBuildings(1, 1) >= wB 
            thetaRightCorner2RIS = atan((q_RIS(2) - S(ii).cornerBuildings(1, 1))/(q_RIS(1) - S(ii).cornerBuildings(1, 2))); % y / x
            thetaLeftCorner2RIS = atan((q_RIS(2) - S(ii).cornerBuildings(2, 1))/(q_RIS(1) - S(ii).cornerBuildings(2, 2)));
            thetaTUAV2RIS = atan2((q_RIS(2) - q_TUAV(2)), (q_RIS(1) - q_TUAV(1)));
            if thetaTUAV2RIS > thetaRightCorner2RIS && thetaTUAV2RIS < thetaLeftCorner2RIS
            %     disp('TUAV is in the cone')
                % In order to know whether the TUAV and RIS has LoS, elevation angles
                % must be compared
                RIS2Building = (q_RIS(1) - S(ii).cornerBuildings(1, 2)) / cos(thetaTUAV2RIS);
                RIS2BuildingElevation = atan((S(ii).height-q_RIS(3))/RIS2Building);
                idxRIS(ii) = RIS2TUAVElevation >= RIS2BuildingElevation;
            end
        end
    
    else
        thetaRightCorner2RIS = atan((q_RIS(2) - S(ii).cornerBuildings(1, 1))/(q_RIS(1) - S(ii).cornerBuildings(1, 2)));
        thetaLeftCorner2RIS = atan((q_RIS(2) - S(ii).cornerBuildings(2, 1))/(q_RIS(1) - S(ii).cornerBuildings(1, 2)));
        thetaTUAV2RIS = atan2((q_RIS(2) - q_TUAV(2)), (q_RIS(1) - q_TUAV(1)));
        if thetaTUAV2RIS > thetaRightCorner2RIS && thetaTUAV2RIS < thetaLeftCorner2RIS
        %     disp('TUAV is in the cone')
            % In order to know whether the TUAV and RIS has LoS, elevation angles
            % must be compared
            RIS2Building = wS / cos(thetaTUAV2RIS);
            RIS2BuildingElevation = atan((S(ii).height-q_RIS(3))/RIS2Building);
            idxRIS(ii) = RIS2TUAVElevation >= RIS2BuildingElevation;
        end
    
    end
end

isRISLOS = all(boolean(idxRIS));
