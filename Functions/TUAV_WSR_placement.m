seed_value = 150; pp = 0; K = 2; N_R = 128; tau = 0.8; precoderIC = 1; 
h_B = 30; % m. It is a 10-storey bulding
q_B = [0, 0, h_B]';
q_TUAV = [0, 0, 0]' + q_B;

initialise_params_differentScenarios;
maximumTetherLength = 100; % m



% linearGrid = -100:5:100;
linearGrid = [-100:10:-10, -8:1:q_RIS(1)];
[X,Y] = meshgrid(linearGrid);
altitudeLevels = 0;
% h_T_U_PL_array = zeros(length(altitudeLevels), K);
h_T_U_PL_array = zeros(length(linearGrid), length(linearGrid), K);
% h_T_R_PL_array = zeros(length(altitudeLevels), 1);
WSR_TUAV_Zposition = zeros([size(X), length(altitudeLevels)]);
for iz = 1:length(altitudeLevels)
    fprintf('iz: %d out of %d\n', iz, length(altitudeLevels))
    altitudeTUAV = altitudeLevels(iz); % altitude from the rooftop
    Z = ones(size(X));
    Z(sqrt(X.^2 + Y.^2 + altitudeTUAV)>maximumTetherLength) = 0;
    for xIndex = 1:length(linearGrid)
        for yIndex = 1:length(linearGrid)
            if Z(xIndex, yIndex) == 1
                q_TUAV = [linearGrid(xIndex), linearGrid(xIndex), altitudeTUAV]' + q_B;
                initialise_params_differentScenarios;
                h_T_U_PL_array(xIndex, yIndex, :) = sum_square_abs(h_T_U_PL);
%                 h_T_R_PL_array(iz, :) = sum(sum_square_abs(G, 2));
                if linearGrid(xIndex) < q_RIS(1)
                    WSR_TUAV_Zposition(xIndex, yIndex, iz) = wsr;
                end
            end
        end
    end
end
disp('finish')

%%
figure;
hold on;
plot3(q_RIS(1), q_RIS(2), q_RIS(3)*0, 'rs')
plot3(q_UEs(1, :), q_UEs(2, :), q_UEs(3, :)*0,'ko')
for k=1:K
    text(q_UEs(1, k)+deltaTxt,q_UEs(2, k)+deltaTxt,sprintf('%d', k))
end
grid on
% legend({'TUAV', 'RIS', 'UEs'})
for iz = 1



    r = rectangle('Position',[-widthBuilding/2 -widthBuilding/2 widthBuilding widthBuilding]);
    rectangle('Position',[-widthBuilding/2+wB -widthBuilding/2 widthBuilding widthBuilding]);
    % surfc(X(WSR_TUAV_Zposition > 0), Y(WSR_TUAV_Zposition > 0), WSR_TUAV_Zposition(WSR_TUAV_Zposition > 0), '.')
    WSR_TUAV_Zposition2  = WSR_TUAV_Zposition(:, :, iz);
    WSR_TUAV_Zposition2(WSR_TUAV_Zposition2 == 0) = NaN;
    surf(Y, X, WSR_TUAV_Zposition2, 'FaceAlpha',0.5');%, 'DisplayName', sprintf('TUAV = %d [m]', altitudeLevels(iz)+h_B))
    text(q_RIS(1), 0, max(WSR_TUAV_Zposition2(:)), sprintf('TUAV = %d [m]', altitudeLevels(iz)))
%     plot3(Y, X, WSR_TUAV_Zposition2, 'b.', 'DisplayName', sprintf('TUAV = %d [m]', altitudeLevels(iz)+h_B))

end

%%

figure; hold on;
for iz = 1%:length(altitudeLevels)
    WSR_TUAV_Zposition2  = WSR_TUAV_Zposition(:, :, iz);
    WSR_TUAV_Zposition2(WSR_TUAV_Zposition2 == 0) = NaN;
    plot(X(Y == 0), WSR_TUAV_Zposition2(X == 0), '--x', 'DisplayName', sprintf('TUAV = %d [m]', altitudeLevels(iz)))

end
leg = legend();
setPlotParams;

%%

figure;
zData = h_T_U_PL_array(:, :, 1);  zData(zData == 0) = NaN;
surf(Y,X, zData); hold on;
zData = h_T_U_PL_array(:, :, 2);  zData(zData == 0) = NaN;
surf(Y,X, zData); 
% surf(Y,X, h_T_U_PL_array(:, :, 2))