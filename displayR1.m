function displayR1(z0, z1, z2, z3, mask, mask1, mask2, mask3)
close(gcf), figure
z0((mask==0)) = nan;
z1(mask1==0) = nan;
z2(mask2==0) = nan;
z3(mask3==0) = nan;

zmin = min(min(z0));
zmax = max(max(z0));

figure(1)
tiledlayout(1, 4);
colormap default;

ax1 = nexttile;
mesh(z0)
% hold on
% x = 1150:1350;
% y = 800*ones(size(x));
% plot3(x, y, z0(800,1150:1350)+0.05, 'r', LineWidth=2);
% x = 1150:1350;
% y = 1200*ones(size(x));
% plot3(x, y, z0(1200,1150:1350)+0.05, 'r', LineWidth=2);
% 
% y = 800:1200;
% x = 1150*ones(size(y));
% plot3(x, y, z0(800:1200,1150)+0.05, 'r', LineWidth=2);
% y = 800:1200;
% x = 1350*ones(size(y));
% plot3(x, y, z0(800:1200,1350)+0.05, 'r', LineWidth=2);
% hold off
zlim([zmin zmax])
xlim tight;
ylim tight;
view(20, 60)
set(gca, 'FontSize', 20);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')
title('(a)', 'FontName','Times New Roman','FontSize',24)

ax2 = nexttile;
mesh(z1)
zlim([zmin zmax])
xlim tight;
ylim tight;
view(20, 60)
set(gca, 'FontSize', 20);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')                                                    
title('(b)', 'FontName','Times New Roman','FontSize',24)

ax3 = nexttile;
mesh(z2)
zlim([zmin zmax])
xlim tight;
ylim tight;
view(20, 60)
set(gca, 'FontSize', 20);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')
title('(c)', 'FontName','Times New Roman','FontSize',24)

ax4 = nexttile;
mesh(z3)
zlim([zmin zmax])
xlim tight;
ylim tight;
view(20, 60)
set(gca, 'FontSize', 20);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')
title('(d)', 'FontName','Times New Roman','FontSize',24)

cb = colorbar;
cb.Label.FontSize = 20;
caxis(ax1,[zmin zmax])
caxis(ax2,[zmin zmax])
caxis(ax3,[zmin zmax])
caxis(ax4,[zmin zmax])
end
