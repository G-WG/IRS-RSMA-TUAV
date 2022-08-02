% cpt = get(groot,'defaultLineLineWidth');
% if cpt < 1
%     set(groot,'defaultLineLineWidth',1.0)
% end
set(groot,'defaultLineLineWidth',2.0)
set(groot,'defaultLineMarkerSize',8);

ax = get(gca);

ax.XLabel.Interpreter = "latex";
ax.XLabel.FontSize = 14;
ax.YLabel.Interpreter = "latex";
ax.YLabel.FontSize = 14;
ax.Title.Interpreter ="latex";
ax.Title.FontSize = 18;
leg.Interpreter = "latex";
leg.FontSize = 14;
leg.Location = "best"; 
grid on;
% axis square;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)
set(gca,'XTickLabelMode','auto')
a = get(gca,'YTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)
set(gca,'XTickLabelMode','auto')