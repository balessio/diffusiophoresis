global parameter_file_upload

parameter_file_unload = matfile(parameter_file_upload);

% nondimensional parameters
Peclet = parameter_file_unload.Peclet;
Dp_Ds= parameter_file_unload.Dp_Ds; % diffusivity of colloids
L_h= parameter_file_unload.L_h; % aspect ratio
n_bulk= parameter_file_unload.n_bulk; % main channel colloid concentration
dx= parameter_file_unload.dx; % step along pore axis
dy= parameter_file_unload.dy; % step along pore width
dt= parameter_file_unload.dt; % time step
t0= parameter_file_unload.t0; % initial time (skipping 0 because its unstable
duration = parameter_file_unload.duration; % end time
z_i = parameter_file_unload.z_i;
D_i = parameter_file_unload.D_i;
PsiWref = parameter_file_unload.PsiWref;
cWref = parameter_file_unload.cWref;
PsiDref = parameter_file_unload.PsiDref;
cDref = parameter_file_unload.cDref;
output_folder = parameter_file_unload.output_folder;
solute_cell_folder = parameter_file_unload.solute_cell_folder;


if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

fid = fopen(join([output_folder,'log.txt'],''), 'wt');
log_string = join([datestr(now, 'mm/dd/yyyy-HH:MM:SS'),'\n'...
    'dx = %f\n','dy = %f\n',...
    'dt = %f\n','start_time = %f\n','duration = %f\n',...
    'n bulk = %f\n','L_h = %f\n','z_i = %f %f %f\n','D_i = %f %f %f\n',...
    'cWref = %f %f %f\n','PsiWref = %f\n',...
    'cDref = %f %f %f\n','PsiDref = %f\n','Peclet = %f\n',...
    'solute cell folder = %f/n'],'');
fprintf(fid, log_string, dx, dy, dt,...
    t0, duration, n_bulk, L_h, z_i, D_i, cWref, PsiWref,...
    cDref, PsiDref, Peclet, solute_cell_folder);
fclose(fid);
save(join([output_folder,'log.mat'],''),...
'dx', 'dy', 'dt',...
't0', 'duration', 'n_bulk', 'L_h', 'z_i', 'D_i',...
'cWref', 'PsiWref','cDref', 'PsiDref','Peclet','solute_cell_folder');

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
    c1_matrices(:,:,i) = solute_cell{i,1}; %#ok<SAGROW>
    c2_matrices(:,:,i) = solute_cell{i,2}; %#ok<SAGROW>
    c3_matrices(:,:,i) = solute_cell{i,3}; %#ok<SAGROW>
end
[x_n, y_n, t_n] = ndgrid(x,y,t);
interpFunc = griddedInterpolant(x_cCell,y_cCell,t_cCell,c1_matrices);
c1_frames = interpFunc(x_n,y_n,t_n);
interpFunc = griddedInterpolant(x_cCell,y_cCell,t_cCell,c2_matrices);
c2_frames = interpFunc(x_n,y_n,t_n);
interpFunc = griddedInterpolant(x_cCell,y_cCell,t_cCell,c3_matrices);
c3_frames = interpFunc(x_n,y_n,t_n);

XYTcolloid_cell=cell(1,3);
XYTcolloid_cell{1,1}=x_c; XYTcolloid_cell{1,2}=y_c; XYTcolloid_cell{1,3}=t;
colloid_cell=cell(length(t),1);

PsiW_array = zeros(length(x),length(y),length(t));
PsiD_array = zeros(length(x),length(y),length(t));

PsiWref = ones(length(x),length(y))*PsiWref;
PsiDref = ones(length(x),length(y))*PsiDref;

tic;
% march in time
for count=1:length(t)
    
    c_array = zeros(length(x),length(y),3);
    c_array(:,:,1) = c1_frames(:,:,count);
    c_array(:,:,2) = c2_frames(:,:,count);
    c_array(:,:,3) = c3_frames(:,:,count);
    
    solve_PsiW = @(PsiW) F_minus_Fref(PsiW,c_array,PsiWref,cWref,z_i);
    PsiW = fsolve(solve_PsiW,PsiWref);
    solve_PsiD = @(PsiD) F_minus_Fref(PsiD,c_array,PsiDref,cDref,z_i);
    PsiD = fsolve(solve_PsiD,PsiDref);
    PsiW_array(:,:,count) = PsiW; PsiD_array(:,:,count) = PsiD;

%     PsiW = -4; PsiD = -3;
%     PsiW_array(:,:,count) = PsiW; PsiD_array(:,:,count) = PsiD;
    
    % velocity evaluate at cell boundaries, concentration in center
    [vp_x,vp_y,vf_x,vf_y]=vp_discrete_multi_ion(c1_frames(:,:,count),...
        c2_frames(:,:,count),c3_frames(:,:,count),D_i,z_i,...
        x,y,dx,dy,Peclet,PsiW,PsiD,L_h);
    
    n0=n;
    vp_x=transpose(vp_x);
    vp_y=transpose(vp_y);
    
    n=fsolve(@(Xval)myfun(Xval,c1_frames(:,:,count),...
        c2_frames(:,:,count),c3_frames(:,:,count),D_i,z_i,...
        x,y,x_c,y_c,dx,dy,t(count),dt,n_bulk,...
        Peclet,PsiW,PsiD,Dp_Ds,L_h,n0),reshape(n0,[length(x_c)*length(y_c),1]));
    n=reshape(n, [length(x_c),length(y_c)]);
    
    [xx_c,yy_c]=meshgrid(x_c,y_c);
    xx_c=transpose(xx_c);
    yy_c=transpose(yy_c);
    
    [xx,yy]=meshgrid(x,y);
    
    colloid_cell{count,1} = n;
    
    VFX = vf_x';
    VFY = vf_y';
    
%     figure(1)
%     drawnow
%     surf(xx_c,yy_c,n)    
%     colorbar;
%     colormap(hot);
%     caxis([0 max(n(:,1/L_h/dy)*1.5)])
%     %caxis([0 3])
%     shading interp
%     axis equal tight;
%     title(strcat('t=',num2str(t(count))))
%     view([0 0 1])
    
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
save(join([output_folder,'potentials.mat'],''), 'PsiW_array', 'PsiD_array');

function [F,f_l,f_r]=myfun(Xval,c1,c2,c3,D_i,z_i,x,y,x_c,y_c,dx,dy,...
    t,dt,n_bulk,Pe,PsiW,PsiD,DpDs,L_h,n0)

n=reshape(Xval,[length(x_c),length(y_c)]);

[vp_x,vp_y,~,~,~,~]=vp_discrete_multi_ion(c1,c2,c3,D_i,z_i,x,y,...
    dx,dy,Pe,PsiW,PsiD,L_h);
% velocity evaluate at cell boundaries

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
f_l(i,j)=  - DpDs*(n(i,j)-n_bulk)/dx/2;

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

%----------------------------------------------
function solve_for_Psi = F_minus_Fref(Psi_,c_,Psi_ref,c_ref,z_i) % solve PsiW with refs
s = size(c_);
c_ref_array = zeros(s);
for ion = 1:length(z_i)
   c_ref_array(:,:,ion) = ones(s(1),s(2))*c_ref(ion);
end
solve_for_Psi = FF(Psi_,c_,z_i)-FF(Psi_ref,c_ref_array,z_i);
end
function Fpsifunc = FF(Psi_,c_,z_i) % solve for Psi with references
s = size(c_);
Fpsifunc = zeros(s(1),s(2));
for i=1:length(z_i)
    Fpsifunc = Fpsifunc + c_(:,:,i).*(exp(-z_i(i).*Psi_) - 1);
end
Fpsifunc = 2*Fpsifunc;
end
%----------------------------------------------
