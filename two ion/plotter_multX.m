function plotter_multX(X,Y,colors,labels,ii,jj)
%plot(X,Y(1,:),'linewidth',3)
box on
dashed_line_index = ii;
dotted_line_index = jj;
for i=1:length(Y)
    %hold on
    if i<=dashed_line_index
    	plot(X{i},Y{i},'linewidth',3,'color',colors(i,:))
    end
    if i>dashed_line_index && i<=dotted_line_index
    	plot(X{i},Y{i},'linewidth',3,'color',colors(i-dashed_line_index,:),'linestyle','--')
    end
    if i>dotted_line_index
    	plot(X{i},Y{i},'linewidth',3,'color',colors(i-dotted_line_index,:),'linestyle',':')
    end
    hold off
    hold on
end
%hold on
l=legend(labels);
legend boxoff
set(l,'fontsize',30,'interpreter','Latex','Location','best')
pbaspect([4 3 1])
set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex')
%set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex','yscale','log')
%set(gca,'linewidth',3,'fontsize',30,'ticklabelinterpreter','latex','xscale','log')
ylim([0 inf])
%xlim([0 1])
ax=gca;
%ax.XTick=[0 0.5 1];
end
