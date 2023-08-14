function displayR2(z0, zMF, zWGIF, zGD, mask)
close(gcf), figure;
z0(mask==0) = nan;
zGD(mask==0) = nan;
zMF(mask==0) = nan;
zWGIF(mask==0) = nan;
zmin = min(min(z0));
zmax = max(max(z0));
figure(1)
tiledlayout(1, 4);
colormap default;

ax1 = nexttile;
mesh(z0)
zlim([zmin zmax])
view(16, 40)
xlim tight
ylim tight
set(gca, 'FontSize', 16);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')
title('(a)', 'FontName','Times New Roman','FontSize',24)

ax2 = nexttile;
mesh(zMF)
zlim([zmin zmax])
view(16, 40)
xlim tight
ylim tight
set(gca, 'FontSize', 16);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')
title('(b)', 'FontName','Times New Roman','FontSize',24)

ax3 = nexttile;
mesh(zWGIF)
zlim([zmin zmax])
view(16, 40)
xlim tight
ylim tight
set(gca, 'FontSize', 16);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')
title('(c)', 'FontName','Times New Roman','FontSize',24)

ax4 = nexttile;
mesh(zGD)
zlim([zmin zmax])
view(16, 40)
xlim tight
ylim tight
set(gca, 'FontSize', 16);
grid off, box on
zlabel('depth (mm)')
ylabel('pixel')
xlabel('pixel')
title('(d)', 'FontName','Times New Roman','FontSize',24)


cb = colorbar;
cb.Label.FontSize = 16;
caxis(ax1,[zmin zmax])
caxis(ax2,[zmin zmax])
caxis(ax3,[zmin zmax])
caxis(ax4,[zmin zmax])
end