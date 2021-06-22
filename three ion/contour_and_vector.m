function contour_and_vector(X,Y,Z,Z_min,Z_max,Vx,Vy)
box on
pcolor(X',Y',Z)    
%colorbar('Ticks',Z_min:(Z_max-Z_min)/2:Z_max+1,'ticklabelinterpreter','latex')
colorbar('Ticks',0:0.55:1.1,'ticklabelinterpreter','latex','fontsize',30)
colormap cool
%caxis([Z_min Z_max])
caxis([0 1.1])
shading interp
axis equal tight;
hold on
%quiver(X',Y',Vx',Vy')
%streamline(X,Y,Vx,Vy,zeros(size(Y(:,1))),Y(:,1))
h = streamslice(X,Y,Vx,Vy,0.5);
set(h, 'color', 'k')
axis([0 1 -0.1 0.1])
%axis([0 1 -0.35 0.35]) %for making large colorbar
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
ax=gca;
ax.XTick=[0 0.25 0.5 0.75 1];
ax.YTick=[-0.1 0 0.1];
view([0 0 1])
end