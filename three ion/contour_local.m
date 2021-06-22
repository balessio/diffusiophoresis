function contour_local(X,Y,Z,Z_min,Z_max)
pcolor(X',Y',Z)    
%colorbar('Ticks',Z_min:(Z_max-Z_min)/2:Z_max+1,'ticklabelinterpreter','latex')
colorbar('Ticks',0:0.95:1.9,'ticklabelinterpreter','latex')
colormap hot
%caxis([Z_min Z_max])
caxis([0 2.0])
shading interp
axis equal tight;
hold on
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
ax=gca;
ax.XTick=[0.02 0.25 0.5 0.75 0.98];
ax.YTick=[-0.095 0 0.095];
view([0 0 1])
end
