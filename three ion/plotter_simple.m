function plotter_simple(X,Y)
%plot(X,Y(1,:),'linewidth',3)
box on
for i=1:size(Y,2)
    %hold on
    plot(X{i},Y{i},'linewidth',3)
    hold off
    hold on
end
end