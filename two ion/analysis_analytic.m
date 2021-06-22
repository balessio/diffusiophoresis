%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colloid_folder = 'data/colloid/test/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colloid_cell_unload = matfile(join([colloid_folder,'colloid_cell.mat'],''));
colloid_cell = colloid_cell_unload.colloid_cell;
XYTcolloid_cell_unload = matfile(join([colloid_folder,'XYTcolloid_cell.mat'],''));
XYTcolloid_cell = XYTcolloid_cell_unload.XYTcolloid_cell;
colloid_log = matfile(join([colloid_folder,'log.mat'],''));

tic;

n_max = 0;
for count=1:length(XYTcolloid_cell{1,3})
    temp_max = max(colloid_cell{count,1}(:,int8(size(colloid_cell{count,1},2)/2)));
    if temp_max > n_max
       n_max = temp_max; 
    end
end
n_max = n_max*1.25;

startyf = [];
startyf(1)=-1/colloid_log.L_h;
for yValIndex = 1:length(XYTcolloid_cell{1,2})
    if rem(yValIndex,2)==0
        startyf = [startyf XYTcolloid_cell{1,2}(yValIndex)]; %#ok<AGROW>
    end
end
%startyf = XYTsolute_cell{1,2};
startxf = zeros(size(startyf));

% march in time
for count=1:length(XYTcolloid_cell{1,3})
    
    % velocity evaluate at cell boundaries, concentration in center
    [vp_x,vp_y,vf_x,vf_y,dcdx,dcdy] = vp(XYTcolloid_cell{1,1},...
        XYTcolloid_cell{1,2},XYTcolloid_cell{1,3}(count),...
        colloid_log.Gp_Ds,colloid_log.Gw_Ds,colloid_log.B,colloid_log.L_h);
    
    [dndx,dndy]=n_flux(colloid_cell{count,1},...
        XYTcolloid_cell{1,1},XYTcolloid_cell{1,2},colloid_log.dx,...
        colloid_log.dy);
    % the X=0 values destroy the scale
    %dndx(1,:)=0; dndy(1,:)=0;
    
    [xx_nCell,yy_nCell]=meshgrid(XYTcolloid_cell{1,1},XYTcolloid_cell{1,2});
    
    f=figure(1);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    drawnow
    subplot_vertical = 5; subplot_horiz = 2;
    ax(1) = subplot(subplot_vertical, subplot_horiz, [1 2]);
    surf(xx_nCell',yy_nCell',colloid_cell{count,1})    
    colorbar('Ticks',0:0.25:n_max+1)
    colormap(ax(1), hot);
    caxis([0 n_max])
    shading interp
    axis equal tight;
    title(strcat('colloid concentration, tau=',num2str(XYTcolloid_cell{1,3}(count))))
    view([0 0 1])
    xlabel('X')
    ylabel('Y')
    
    ax(2) = subplot(subplot_vertical, subplot_horiz, [3 4 5 6]);
    quiver(xx_nCell',yy_nCell',vf_x',vf_y')
    %streamline(xx_nCell,yy_nCell,vf_x,vf_y,startxf,startyf)
    title('fluid velocity vectors [Vfx, Vfy]')
    axis([0 1 -1/colloid_log.L_h 1/colloid_log.L_h])
    xlabel('X')
    ylabel('Y')
    
%     ax(4) = subplot(subplot_vertical, subplot_horiz, [7 8]);
%     quiver(xx_nCell',yy_nCell',dndx,dndy)
%     %startxn = XYTcolloid_cell{1,1}; startyn = zeros(size(startxn));
%     %streamline(xx_nCell,yy_nCell,dndx',dndy',startxn,startyn)
%     title('flux vectors [dn/dx, dn/dy], X=0 values eliminated')
%     axis([0 1 -1/colloid_log.L_h 1/colloid_log.L_h])
%     xlabel('X')
%     ylabel('Y')
    
    ax(3) = subplot(subplot_vertical, subplot_horiz, [7 9]);
    plot(XYTcolloid_cell{1,1},mean(colloid_cell{count,1},2),'-k','linewidth',1.5)
    axis([0 1 0 n_max])
    title('Y-averaged colloid concentration vs X')
    xlabel('X')
    ylabel('< n >')
    grid on 
    
%     ax(4) = subplot(subplot_vertical, subplot_horiz, [9 11]);
%     plot(XYTsolute_cell{1,1},mean(vp_x,1),'-k','linewidth',1.5)
%     axis([0 1 0 n_max])
%     title('Y-averaged colloid x-velocity vs X')
%     xlabel('X')
%     ylabel('< v_{px} >')
%     grid on
    
    ax(4) = subplot(subplot_vertical, subplot_horiz, [8 10]);
    plot(XYTcolloid_cell{1,1},...
        colloid_cell{count,1}(:,int8(size(colloid_cell{count,1},2)/2)),'-k','linewidth',1.5)
    axis([0 1 0 n_max])
    title('colloid concentration vs X at Y=0')
    xlabel('X')
    ylabel('n')
    grid on
    
    Fvid(count) = getframe(gcf); %#ok<SAGROW>
    clf(f);
    
    int_x = trapz(XYTcolloid_cell{1,1},colloid_cell{count,1},1);
    int_xy = trapz(XYTcolloid_cell{1,2},int_x);
    num_particles(count) = int_xy; %#ok<SAGROW>

end
toc;

if ~exist(join([colloid_folder,'analysis'],''), 'dir')
    mkdir(join([colloid_folder,'analysis'],''))
end
writerObj = VideoWriter(join([colloid_folder,'analysis/movie.avi'],''));
writerObj.FrameRate = 16;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(Fvid)
    % convert the image to a frame
    frame = Fvid(i);
    writeVideo(writerObj, frame);
end
close(writerObj);

figure(2)
plot(XYTcolloid_cell{1,3},num_particles)    
xlabel('tau')
ylabel('normalized number of particles')
grid on
saveas(gcf,join([colloid_folder,'analysis/num_particles.png'],''));

function [dndx,dndy]=n_flux(n,x,y,dx,dy)

dndx = [(diff(n)/dx); zeros(1,length(y))];
dndy = [(diff(n')/dy)', zeros(length(x),1)];

end
