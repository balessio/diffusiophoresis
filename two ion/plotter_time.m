function plotter(X,Y,colors,labels)
%plot(X,Y(1,:),'linewidth',3)
box on
for i=1:size(Y,1)
    %hold on
    plot(X,Y(i,:),'linewidth',3,'color',colors(i,:))
    hold off
    hold on
end
%hold on
l=legend(labels);
legend boxoff
set(l,'fontsize',20,'interpreter','Latex','Location','Best')
pbaspect([4 3 1])
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
%set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex','yscale','log')
%set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex','xscale','log')
ax=gca;
t_min = round(min(X),1); t_max = round(max(X),1);
ax.XTick=[t_min (t_max-t_min)/2 t_max];
%xlim([0 2])
end