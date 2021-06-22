global parameter_file_upload
parameter_file_unload = matfile(parameter_file_upload);

Peclet = parameter_file_unload.Peclet;
PsiD = parameter_file_unload.PsiD;
n_bulk = parameter_file_unload.n_bulk;
B = parameter_file_unload.B;
L_h = parameter_file_unload.L_h;
PsiW = parameter_file_unload.PsiW;
D_pos_D_Na = parameter_file_unload.D_pos_D_Na;
D_neg_D_Na = parameter_file_unload.D_neg_D_Na;
Dp_D_Na = parameter_file_unload.Dp_D_Na;
t0 = parameter_file_unload.t0;
duration = parameter_file_unload.duration;
dx = parameter_file_unload.dx;
dy = parameter_file_unload.dy;
dt = parameter_file_unload.dt;
solute_cell_folder = parameter_file_unload.solute_folder;
output_folder = parameter_file_unload.output_folder;

Gw_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiW+PsiW^2/8);
Gp_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiD+PsiD^2/8);

global c_analytic
c_analytic = false;

Dp_Ds = Dp_D_Na;

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

save(join([output_folder,'log.mat'],''),...
'solute_cell_folder','Gp_Ds','Gw_Ds','Dp_Ds','L_h','B',...
'n_bulk','dx','dy','dt','t0','duration');

% location of cell boundaries
x=0:dx:1;
y=-1/L_h:dy:1/L_h;
t=t0:dt:duration;

% location of cell centers
x_c=0.5*(x(2:end)+x(1:end-1));
y_c=0.5*(y(2:end)+y(1:end-1));

% initialize n
n=ones(length(x_c),length(y_c));
f_l=n.*0;
f_r=n.*0;
f_b=n.*0;
f_t=n.*0;

solute_cell_unload = matfile(join([solute_cell_folder,'solute_cell.mat'],''));
solute_cell = solute_cell_unload.solute_cell;
XYTsolute_cell_unload = matfile(join([solute_cell_folder,'XYTsolute_cell.mat'],''));
XYTsolute_cell = XYTsolute_cell_unload.XYTsolute_cell;
[x_cCell, y_cCell, t_cCell] = ndgrid(XYTsolute_cell{1,1},...
    XYTsolute_cell{1,2},XYTsolute_cell{1,3});
c_matrices = [];
for i=1:length(solute_cell)
    c_matrices(:,:,i) = solute_cell{i,1}; %#ok<SAGROW>
end
interpFunc = griddedInterpolant(x_cCell,y_cCell,t_cCell,c_matrices);
[x_n, y_n, t_n] = ndgrid(x,y,t);
c_frames = interpFunc(x_n,y_n,t_n);

XYTcolloid_cell=cell(1,3);
XYTcolloid_cell{1,1}=x_c; XYTcolloid_cell{1,2}=y_c; XYTcolloid_cell{1,3}=t;
colloid_cell=cell(length(t),1);

tic;
% march in time
for count=1:length(t)
    
    % velocity evaluate at cell boundaries, concentration in center
    if c_analytic
        [vp_x,vp_y,~,~,~,~] = vp(x,y,t(count),Gp_Ds,Gw_Ds,B,L_h);
    else
        [vp_x,vp_y,~,~,~,~,~,~] = vp_discrete(c_frames(:,:,count),...
            x,y,dx,dy,Gp_Ds,Gw_Ds,L_h);
    end
    
    n0=n;
    vp_x=transpose(vp_x);
    vp_y=transpose(vp_y);
    
    n=fsolve(@(Xval)myfun(Xval,c_frames(:,:,count),x,y,x_c,y_c,dx,dy,...
        t(count),dt,B,n_bulk,...
        Gp_Ds,Gw_Ds,Dp_D_Na,L_h,n0),reshape(n0,[length(x_c)*length(y_c),1]));
    n=reshape(n, [length(x_c),length(y_c)]);
    
    [xx_c,yy_c]=meshgrid(x_c,y_c);
    xx_c=transpose(xx_c);
    yy_c=transpose(yy_c);
    
    [xx,yy]=meshgrid(x,y);
    
    colloid_cell{count,1} = n;
    
    figure(1)
    drawnow
    ax(1) = subplot(2, 1, 1);
    surf(xx_c,yy_c,n)    
    colorbar;
    colormap(hot);
    caxis([0 max(n(:,1/L_h/dy)*1.5)])
    %caxis([0 3])
    shading interp
    axis equal tight;
    title(strcat('t=',num2str(t(count))))
    view([0 0 1])
    
    ax(2) = subplot(2, 1, 2);
    %plot(x_c,mean(n,2),'-k','linewidth',1.5)
    plot(x_c,n(:,int8(size(n,2)/2)),'-k','linewidth',1.5)
    axis([0 1 0 max(n(:,1/L_h/dy)*1.5)])
    %title('Y-averaged colloid concentration vs X')
    title('center line colloid concentration vs X')
    xlabel('X')
    %ylabel('< n >')
    ylabel('n(Y=0)')
    grid on

end
toc;

save(join([output_folder,'XYTcolloid_cell.mat'],''), 'XYTcolloid_cell');
save(join([output_folder,'colloid_cell.mat'],''), 'colloid_cell');

function [F,f_l,f_r]=myfun(Xval,c,x,y,x_c,y_c,dx,dy,t,dt,B,n_bulk,...
    GpDs,GwDs,DpDs,L_h,n0)
global c_analytic

n=reshape(Xval,[length(x_c),length(y_c)]);

% velocity evaluate at cell boundaries
if c_analytic
    [vp_x,vp_y,~,~,~,~]=vp(x,y,t,GpDs,GwDs,B,L_h);
else
    [vp_x,vp_y,~,~,~,~]=vp_discrete(c,x,y,dx,dy,GpDs,GwDs,L_h);
end

vp_x=transpose(vp_x);
vp_y=transpose(vp_y);

%%%%% HORIZONTAL FLUXES
j=1:length(y_c);

% left flux
i=2:length(x_c);
vp_l = (vp_x(i,j)+vp_x(i,j+1))*0.5;
n_l = (vp_l>0).*n(i-1,j) + (vp_l <= 0).*n(i,j);
f_l(i,j) = vp_l.*n_l - DpDs*(n(i,j)-n(i-1,j))/dx;
i=1;
f_l(i,j)= - DpDs*(n(i,j)-n_bulk)/dx/2;

% right flux
i=1:length(x_c)-1;
vp_r = (vp_x(i+1,j)+vp_x(i+1,j+1))*0.5;
n_r = (vp_r>0).*n(i,j) + (vp_r <= 0).*n(i+1,j);
f_r(i,j) = vp_r.*n_r - DpDs*(n(i+1,j)-n(i,j))/dx;
i=length(x_c);
f_r(i,j)=0;

%%%%% VERTICAL FLUXES

i=1:length(x_c);
%bottom flux
j=2:length(y_c);
vp_b = (vp_y(i,j)+vp_y(i+1,j))*0.5;
n_b = (vp_b>0).*n(i,j-1) + (vp_b <= 0).*n(i,j);
f_b(i,j) = vp_b.*n_b - DpDs*(n(i,j)-n(i,j-1))/dy;
j=1;
f_b(i,j)=0;

%top flux
j=1:length(y_c)-1;
vp_t = (vp_y(i,j+1)+vp_y(i+1,j+1))*0.5;
n_t = (vp_t>0).*n(i,j) + (vp_t <= 0).*n(i,j+1);
f_t(i,j) = vp_t.*n_t - DpDs*(n(i,j+1)-n(i,j))/dy;

j=length(y_c);
f_t(i,j)=0;

F = reshape( (n-n0 - dt*((f_l-f_r)/dx + (f_b-f_t)/dy)),...
    [length(x_c)*length(y_c),1] );

end
