maxPosTUAV2ndScenery = linearGrid(maxPoints_2s(1:2, 1));
maxPosTUAV1stScenery = linearGrid(maxPoints(1:2, 1));

Pt_dB_list = [0:5:50];
WSR_Pt = zeros(length(Pt_dB_list), 1);
Rate = zeros(length(Pt_dB_list), 4); % rate_kp, rate_c
for iPt = 1:length(Pt_dB_list)
    Pt_dB = Pt_dB_list(iPt);
    %% Scenary
    q_TUAV = [maxPosTUAV1stScenery, h_B]';
    initialise_params_differentScenarios_v2;
    
    % q_TUAV = maxPosTUAV2ndScenery;
    % initialise_params_differentScenarios2nd_v2;

    %% WSR
    
    WSR_Pt(iPt) = wsr;
    Rate(iPt, :) = [rate_kp', rate_c'];
end

%%

% figure;


