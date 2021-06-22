function plot_one_line_time(X,Y)
plot(X,Y,'-b','linewidth',3)
hold on
pbaspect([4 3 1])
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
ax=gca;
t_min = round(min(X),1); t_max = round(max(X),1);
ax.XTick=[t_min (t_max-t_min)/2 t_max];
end