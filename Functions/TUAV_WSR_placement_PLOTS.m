% %%
% legbool =1;
% deltaTxt = 0.5;
% figure;
% % plot3(q_TUAV(1), q_TUAV(2), q_TUAV(3), '*'); 
% subplot(1, 2, 1)
% hold on;
% q_UEs = zeros(3, K);
% q_UEs(:, 1) = [q_RIS(1)-3, q_RIS(2)+3, 1.5];
% q_UEs(:, 2) = [q_RIS(1)-wS+1, q_RIS(2)-1, 1.5];
% plot3(q_RIS(1), q_RIS(2), q_RIS(3), 'rs')
% plot3(q_UEs(1, :), q_UEs(2, :), q_UEs(3, :),'ko')
% for k=1:K
%     text(q_UEs(1, k)+deltaTxt,q_UEs(2, k)+deltaTxt,sprintf('%d', k))
% end
% grid on
% leg = legend({'RIS', 'UEs'})
% 
% r = rectangle('Position',[-widthBuilding/2 -widthBuilding/2 widthBuilding widthBuilding]);
% rectangle('Position',[-widthBuilding/2+wB -widthBuilding/2 widthBuilding widthBuilding]);
% xlabel('[m]')
% ylabel('[m]')
% xlim([-10, 30])
% ylim([-20, 20])
% title('Scenery-1')
% setPlotParams;
% 
% subplot(1, 2, 2); hold on;
% q_UEs = zeros(3, K);
% q_UEs(:, 1) = [q_RIS(1)-3, q_RIS(2)+3, 1.5];
% q_UEs(:, 2) = [-wS/2-1, q_RIS(2)-1, 1.5];
% plot3(q_RIS(1), q_RIS(2), q_RIS(3), 'rs')
% plot3(q_UEs(1, :), q_UEs(2, :), q_UEs(3, :),'ko')
% for k=1:K
%     text(q_UEs(1, k)+deltaTxt,q_UEs(2, k)+deltaTxt,sprintf('%d', k))
% end
% grid on
% leg = legend({'RIS', 'UEs'})
% 
% r = rectangle('Position',[-widthBuilding/2 -widthBuilding/2 widthBuilding widthBuilding]);
% rectangle('Position',[-widthBuilding/2+wB -widthBuilding/2 widthBuilding widthBuilding]);
% xlabel('[m]')
% ylabel('[m]')
% xlim([-10, 30])
% ylim([-20, 20])
% title('Scenery-2')
% setPlotParams;
%% Max WSR per altitude Layer
maxPoints = zeros(4, length(altitudeLevels));
for iz=1:length(altitudeLevels)
    WSR_TUAV_Zposition2  = WSR_TUAV_Zposition(:, :, iz);
    [a, b] = max(WSR_TUAV_Zposition2(:));
    [I1,I2,I3,I4] = ind2sub(size(WSR_TUAV_Zposition2),b);
    maxPoints(:, iz) = [I1, I2, altitudeLevels(iz), a]; % [x, y, z, wsr]
end

% maxPoints_2s = zeros(4, length(altitudeLevels));
% for iz=1:length(altitudeLevels)
%     WSR_TUAV_Zposition2  = WSR_TUAV_Zposition_2nd(:, :, iz);
%     [a, b] = max(WSR_TUAV_Zposition2(:));
%     [I1,I2,I3,I4] = ind2sub(size(WSR_TUAV_Zposition2),b);
%     maxPoints_2s(:, iz) = [I1, I2, altitudeLevels(iz), a]; % [x, y, z, wsr]
% end

%%
figure; 
hold on;
o1 = plot(maxPoints(3, :)+h_B, maxPoints(4, :)*B/1e6, '--x');
% o2 = plot(maxPoints_2s(3, :)+h_B, maxPoints_2s(4, :)*B/1e6, '--o');
o3 = plot([h_B, h_B], [0, 550], 'k--');
xlabel('TUAV Altitude [m]')
ylabel('WSR [Mbps]')
ylim([min(maxPoints(4, :)*B/1e6), max(maxPoints(4, :)*B/1e6)])
% leg = legend('Scenery-1', 'Scenery-2', 'Building height');
legbool = 0;
% ylim([250, 550])
setPlotParams;
axis square;

%% WSR 3D

d = 3;
figure;
hold on;
o1 = plot3(q_RIS(1), q_RIS(2), q_RIS(3)*0, 'rs');
o2 = plot3(q_UEs(1, :), q_UEs(2, :), q_UEs(3, :)*0,'ko');
% for k=1:K
%     text(q_UEs(1, k)+deltaTxt,q_UEs(2, k)+deltaTxt,sprintf('%d', k))
% end
% grid on
leg = legend([o1, o2], {'RIS', 'UEs'});

iz = 2;
WSR_TUAV_Zposition2  = WSR_TUAV_Zposition(1:d:end, 1:d:end, iz);
WSR_TUAV_Zposition2(WSR_TUAV_Zposition2 == 0) = NaN;
surf(Y(1:d:end, 1:d:end), X(1:d:end, 1:d:end), WSR_TUAV_Zposition2*B/1e6, 'FaceAlpha',0.5);%, 'DisplayName', sprintf('TUAV = %d [m]', altitudeLevels(iz)+h_B))
colorbar;

text(q_RIS(1), 0, max(WSR_TUAV_Zposition2(:)), sprintf('TUAV = %d [m]', altitudeLevels(iz)+h_B))

for xi = [-2, xiList]
    for yi = yiList
        rectangle('Position',[-widthBuilding/2-xi*widthBuilding ...
            -widthBuilding/2-stretchingFactor-yi*(widthBuilding+stretchingFactor) ...
            widthBuilding widthBuilding+stretchingFactor*2]);
        xc = -widthBuilding/2-xi*widthBuilding;
        yc = -widthBuilding/2-stretchingFactor-yi*(widthBuilding+stretchingFactor);
%         if xi > -2
%             text(xc+deltaTxt,yc+deltaTxt,sprintf('%d', cpt))
%             cpt = cpt+1;
%         end
    end
end
xlim([-100, 100])
ylim([-100, 100])

xlabel('[m]')
ylabel('[m]')
zlabel('WSR [Mbps]')
legbool = 0;
setPlotParams;
axis square;


%%
% namePlot = 'surfWSR_TUAV-10m_Cityscenery_RSMA_IC_AMBF_MRT';
% namePlot = 'surfWSR_TUAV-10m_Cityscenery_RSMA_20JointIter_5WMMSEmaxIter';
% filePath = '.\IRS-RSMA-TUAV\plots\WSR_TUAV-ALtitude\v2\';
% saveFiguresMultipleFormats(namePlot, filePath)
