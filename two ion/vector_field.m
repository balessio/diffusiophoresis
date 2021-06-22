function vector_field(X,Y,Vx,Vy)
shading interp
axis equal tight;
box on
hold on
%quiver(X',Y',Vx',Vy')
%streamline(X,Y,Vx,Vy,zeros(size(Y(:,1))),Y(:,1))
h = streamslice(X,Y,Vx,Vy,0.5);
set(h, 'color', 'k')
set(h, 'linewidth',3)
%axis([0 1 -0.1 0.1])
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
%ax=gca;
%ax.XTick=[0 0.25 0.5 0.75 1];
%ax.YTick=[-0.1 0 0.1];
view([0 0 1])
end