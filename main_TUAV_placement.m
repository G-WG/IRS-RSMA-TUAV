seed_value = 150; pp = 0; K = 2; N_R = 128; tau = 0.8; precoderIC = 1; 
h_B = 30; % m. It is a 10-storey bulding
q_B = [0, 0, h_B]';
q_TUAV = [0, 0, 0]' + q_B;
initialise_params_differentScenarios2nd_v4;
maximumTetherLength = 100; % m

linearGrid = -maximumTetherLength:1:maximumTetherLength;
% linearGrid = [-100:10:-20, -18:0.5:18, 20:10:100];
[X,Y] = meshgrid(linearGrid);
altitudeLevels = [0:10:100];
WSR_TUAV_Zposition = zeros([size(X), length(altitudeLevels)]);

for iz = 1:length(altitudeLevels)
    
    fprintf('iz: %d out of %d\n', iz, length(altitudeLevels))
    altitudeTUAV = altitudeLevels(iz); % altitude from the rooftop
    
    Z = ones(size(X));
    Z(sqrt(X.^2 + Y.^2 + altitudeTUAV.^2)>maximumTetherLength) = 0;
    
    for xIndex = 1:length(linearGrid)
        for yIndex = 1:length(linearGrid)
            if Z(xIndex, yIndex) == 1
                q_TUAV = [linearGrid(xIndex), linearGrid(yIndex), altitudeTUAV]' + q_B;
                if linearGrid(xIndex) <= q_RIS(1)
                    tic;
                    main;
                    timeElapsed = toc;
                    fprintf('(%d, %d) of %d, timeElapsed %.2f [s]\n', xIndex, yIndex, length(linearGrid), timeElapsed);
                    WSR_TUAV_Zposition(xIndex, yIndex, iz) = wsr;
                end
            end
        end
    end
end
disp('finish')
if precoderIC == 1
    name = sprintf('Hovering_JointOptimisation_RSMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_AMBF_RZF', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R);
else
    name = sprintf('Hovering_JointOptimisation_RSMA-WMMSE_PS-SCAref_%dRth%d_NR_%d_SVD_MRT', floor(Rth(1)), round(mod(Rth(1), 1)*10), N_R);
end
name_matfile = sprintf('logs/JointOptimisation/log_seed_%d_%s.mat', seed_value, name);
save(name_matfile)