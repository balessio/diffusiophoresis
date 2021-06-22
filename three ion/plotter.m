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
set(l,'fontsize',30,'interpreter','latex','Location','Best')
pbaspect([4 3 1])
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
ax=gca;
ax.XTick=[0 0.5 1];
end