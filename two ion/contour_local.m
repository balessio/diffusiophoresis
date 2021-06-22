function contour_local(X,Y,Z,Z_min,Z_max)
pcolor(X',Y',Z)    
colorbar('Ticks',Z_min:(Z_max-Z_min)/2:Z_max+1,'ticklabelinterpreter','latex')
%colorbar('Ticks',0:0.95:1.9,'ticklabelinterpreter','latex')
%colorbar('Ticks',0:2.1:4.2,'ticklabelinterpreter','latex')
%caxis([0 4.2])
%colorbar('Ticks',-2:2:2,'ticklabelinterpreter','latex')
%caxis([-2 2])
%colorbar('Ticks',0:7.5:15,'ticklabelinterpreter','latex')
%caxis([0 15])
colormap hot
%caxis([Z_min Z_max])
%caxis([0 2.0])
shading interp
axis equal tight;
hold on
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
ax=gca;
ax.XTick=[0.02 0.25 0.5 0.75 0.98];
%ax.YTick=[-0.095 0 0.095];
ax.YTick=[-0.046 0 0.046];
view([0 0 1])
end
