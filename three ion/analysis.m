%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global colloid_folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colloid_cell_unload = matfile(join([colloid_folder,'colloid_cell.mat'],''));
colloid_cell = colloid_cell_unload.colloid_cell;
XYTcolloid_cell_unload = matfile(join([colloid_folder,'XYTcolloid_cell.mat'],''));
XYTcolloid_cell = XYTcolloid_cell_unload.XYTcolloid_cell;

potentials_unload = matfile(join([colloid_folder,'potentials.mat'],''));
PsiW_array = potentials_unload.PsiW_array;
PsiD_array = potentials_unload.PsiD_array;

colloid_log = matfile(join([colloid_folder,'log.mat'],''));
solute_log = matfile(join([colloid_log.solute_cell_folder,'log.mat'],''));
solute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'solute_cell.mat'],''));
solute_cell = solute_cell_unload.solute_cell;
XYTsolute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'XYTsolute_cell.mat'],''));
XYTsolute_cell = XYTsolute_cell_unload.XYTsolute_cell;

[x_PsiInterpBefore, y_PsiInterpBefore, t_PsiInterpBefore] = ndgrid(0:colloid_log.dx:1,...
    -1/colloid_log.L_h:colloid_log.dy:1/colloid_log.L_h,...
    XYTcolloid_cell{1,3});
[x_PsiInterpAfter, y_PsiInterpAfter, t_PsiInterpAfter] = ndgrid(...
    XYTsolute_cell{1,1},XYTsolute_cell{1,2},XYTcolloid_cell{1,3});

interpFunc = griddedInterpolant(x_PsiInterpBefore,...
    y_PsiInterpBefore,t_PsiInterpBefore,PsiW_array);
PsiW_interpolated = interpFunc(x_PsiInterpAfter, y_PsiInterpAfter, t_PsiInterpAfter);

interpFunc = griddedInterpolant(x_PsiInterpBefore,...
    y_PsiInterpBefore,t_PsiInterpBefore,PsiD_array);
PsiD_interpolated = interpFunc(x_PsiInterpAfter, y_PsiInterpAfter, t_PsiInterpAfter);

c1_matrices = [];
for i=1:length(solute_cell)
    c1_matrices(:,:,i) = solute_cell{i,1}; %#ok<SAGROW>
end

c2_matrices = [];
for i=1:length(solute_cell)
    c2_matrices(:,:,i) = solute_cell{i,2}; %#ok<SAGROW>
end

c3_matrices = [];
for i=1:length(solute_cell)
    c3_matrices(:,:,i) = solute_cell{i,3}; %#ok<SAGROW>
end

[x_soluteInterpBefore, y_soluteInterpBefore, t_soluteInterpBefore] = ndgrid(...
    XYTsolute_cell{1,1},XYTsolute_cell{1,2},XYTsolute_cell{1,3});
[x_soluteInterpAftere, y_soluteInterpAfter, t_soluteInterpAfter] = ndgrid(...
    XYTsolute_cell{1,1},XYTsolute_cell{1,2},XYTcolloid_cell{1,3});

interpFunc = griddedInterpolant(x_soluteInterpBefore,...
    y_soluteInterpBefore,t_soluteInterpBefore,c1_matrices);
c1_matrices = interpFunc(x_soluteInterpAftere, y_soluteInterpAfter, t_soluteInterpAfter);
interpFunc = griddedInterpolant(x_soluteInterpBefore,...
    y_soluteInterpBefore,t_soluteInterpBefore,c2_matrices);
c2_matrices = interpFunc(x_soluteInterpAftere, y_soluteInterpAfter, t_soluteInterpAfter);
interpFunc = griddedInterpolant(x_soluteInterpBefore,...
    y_soluteInterpBefore,t_soluteInterpBefore,c3_matrices);
c3_matrices = interpFunc(x_soluteInterpAftere, y_soluteInterpAfter, t_soluteInterpAfter);

tic;

n_max = 0;
for count=int64(length(XYTcolloid_cell{1,3})/2):length(XYTcolloid_cell{1,3})
    temp_max = max(colloid_cell{count,1}(:,int8(size(colloid_cell{count,1},2)/2)));
    if temp_max > n_max
       n_max = temp_max; 
    end
end
n_max = n_max*1.25;

startyf = [];
startyf(1)=-1/solute_log.L_h;
for yValIndex = 1:length(XYTsolute_cell{1,2})
    if rem(yValIndex,2)==0
        startyf = [startyf XYTsolute_cell{1,2}(yValIndex)]; %#ok<AGROW>
    end
end
%startyf = XYTsolute_cell{1,2};
startxf = zeros(size(startyf));

c1_max = max(max(max(c1_matrices(:,:,:))));
c2_max = max(max(max(c2_matrices(:,:,:))));
c3_max = max(max(max(c3_matrices(:,:,:))));
c1_min = min(min(min(c1_matrices(:,:,:))));
c2_min = min(min(min(c2_matrices(:,:,:))));
c3_min = min(min(min(c3_matrices(:,:,:))));

% march in time
for count=1:length(XYTcolloid_cell{1,3})
    
    % velocity evaluate at cell boundaries, concentration in center
    [vp_x,vp_y,vf_x,vf_y,vDp_x,vDp_y]=vp_discrete_multi_ion(c1_matrices(:,:,count),...
        c2_matrices(:,:,count),c3_matrices(:,:,count),...
        colloid_log.D_i,colloid_log.z_i,XYTsolute_cell{1,1},...
        XYTsolute_cell{1,2},colloid_log.dx,colloid_log.dy,...
        colloid_log.Peclet,PsiW_interpolated(:,:,count),...
        PsiD_interpolated(:,:,count),solute_log.L_h);
    
    [dndx,dndy]=n_flux(colloid_cell{count,1},...
        XYTcolloid_cell{1,1},XYTcolloid_cell{1,2},colloid_log.dx,...
        colloid_log.dy);
    % the X=0 values destroy the scale
    %dndx(1,:)=0; dndy(1,:)=0;
    
    [xx_cCell,yy_cCell]=meshgrid(XYTsolute_cell{1,1},XYTsolute_cell{1,2});
    [xx_nCell,yy_nCell]=meshgrid(XYTcolloid_cell{1,1},XYTcolloid_cell{1,2});
    
    f=figure(1);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    drawnow
    subplot_vertical = 8; subplot_horiz = 2;
    ax(1) = subplot(subplot_vertical, subplot_horiz, [1 2]);
    surf(xx_nCell',yy_nCell',colloid_cell{count,1})    
    colorbar('Ticks',0:n_max/4:n_max+1)
    colormap(ax(1), hot);
    caxis([0 n_max])
    shading interp
    axis equal tight;
    title(strcat('colloid concentration, tau=',num2str(XYTcolloid_cell{1,3}(count))))
    view([0 0 1])
    xlabel('X')
    ylabel('Y')
    
    ax(2) = subplot(subplot_vertical, subplot_horiz, [3 4]);
    surf(xx_cCell',yy_cCell',c1_matrices(:,:,count))    
    colorbar('Ticks',c1_min:(c1_max-c1_min)/4:c1_max+1)
    colormap(ax(2), cool);
    caxis([c1_min c1_max])
    shading interp
    axis equal tight;
    title(strcat('solute concentration (mM), tau=',num2str(XYTcolloid_cell{1,3}(count))))
    view([0 0 1])
    xlabel('X')
    ylabel('Y')
    
    ax(3) = subplot(subplot_vertical, subplot_horiz, [5 6]);
    surf(xx_cCell',yy_cCell',c2_matrices(:,:,count))    
    colorbar('Ticks',c2_min:(c2_max-c2_min)/4:c2_max+1)
    colormap(ax(3), cool);
    caxis([c2_min c2_max])
    shading interp
    axis equal tight;
    title(strcat('solute concentration (mM), tau=',num2str(XYTcolloid_cell{1,3}(count))))
    view([0 0 1])
    xlabel('X')
    ylabel('Y')
    
    ax(4) = subplot(subplot_vertical, subplot_horiz, [7 8]);
    surf(xx_cCell',yy_cCell',c3_matrices(:,:,count))    
    colorbar('Ticks',c3_min:(c3_max-c3_min)/4:c3_max+1)
    colormap(ax(4), cool);
    caxis([0 c3_max])
    shading interp
    axis equal tight;
    title(strcat('solute concentration (mM), tau=',num2str(XYTcolloid_cell{1,3}(count))))
    view([0 0 1])
    xlabel('X')
    ylabel('Y')
    
    ax(5) = subplot(subplot_vertical, subplot_horiz, [9 10]);
    quiver(xx_cCell',yy_cCell',vf_x',vf_y')
    %streamline(xx_cCell,yy_cCell,vf_x,vf_y,startxf,startyf)
    title('fluid velocity vectors [Vfx, Vfy]')
    axis([0 1 -1/colloid_log.L_h 1/colloid_log.L_h])
    xlabel('X')
    ylabel('Y')
    
    ax(6) = subplot(subplot_vertical, subplot_horiz, [11 12]);
    quiver(transpose(xx_cCell(1:3:end)),transpose(yy_cCell(1:3:end)),...
        transpose(vDp_x(1:3:end)),transpose(vDp_y(1:3:end)))
    %streamline(xx_cCell,yy_cCell,vf_x,vf_y,startxf,startyf)
    title('diffusiophoretic velocity vectors [Vdpx, Vdpy]')
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
    
    ax(7) = subplot(subplot_vertical, subplot_horiz, [13 15]);
    plot(XYTcolloid_cell{1,1},mean(colloid_cell{count,1},2),'-k','linewidth',1.5)
    axis([0 1 0 n_max])
    title('Y-averaged colloid concentration vs X')
    xlabel('X')
    ylabel('< n >')
    grid on
    
%     ax(6) = subplot(subplot_vertical, subplot_horiz, [13 15]);
%     plot(XYTsolute_cell{1,1},mean(vp_x,1),'-k','linewidth',1.5)
%     axis([0 1 0 n_max])
%     title('Y-averaged colloid x-velocity vs X')
%     xlabel('X')
%     ylabel('< v_{px} >')
%     grid on
    
    ax(8) = subplot(subplot_vertical, subplot_horiz, [14 16]);
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
