function plot_one_line(X,Y)
plot(X,Y,'-b','linewidth',3)
hold on
pbaspect([4 3 1])
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
%set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex','xscale','log')
ax=gca;
%ax.XTick=[0 0.5 1];
end