%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global colloid_folder which_time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colloid_cell_unload = matfile(join([colloid_folder,'colloid_cell.mat'],''));
colloid_cell = colloid_cell_unload.colloid_cell;
XYTcolloid_cell_unload = matfile(join([colloid_folder,'XYTcolloid_cell.mat'],''));
XYTcolloid_cell = XYTcolloid_cell_unload.XYTcolloid_cell;
colloid_log = matfile(join([colloid_folder,'log.mat'],''));
solute_log = matfile(join([colloid_log.solute_cell_folder,'log.mat'],''));
solute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'solute_cell.mat'],''));
solute_cell = solute_cell_unload.solute_cell;
XYTsolute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'XYTsolute_cell.mat'],''));
XYTsolute_cell = XYTsolute_cell_unload.XYTsolute_cell;

time_integrate_range = ind_closest(which_time,XYTcolloid_cell{1,3});

c_matrices = [];
for i=1:length(solute_cell)
    c_matrices(:,:,i) = solute_cell{i,1}; %#ok<SAGROW>
end
[x_interpGrid, y_interpGrid, t_interpGrid] = ndgrid(XYTsolute_cell{1,1},...
    XYTsolute_cell{1,2},XYTsolute_cell{1,3});
interpFunc = griddedInterpolant(x_interpGrid,...
    y_interpGrid,t_interpGrid,c_matrices);
[x_interp, y_interp, t_interp] = ndgrid(...
    XYTsolute_cell{1,1},XYTsolute_cell{1,2},XYTcolloid_cell{1,3});
c_interpolated = interpFunc(x_interp, y_interp, t_interp);

startyf = [];
startyf(1)=-1/solute_log.L_h;
for yValIndex = 1:length(XYTsolute_cell{1,2})
    if rem(yValIndex,2)==0
        startyf = [startyf XYTsolute_cell{1,2}(yValIndex)]; %#ok<AGROW>
    end
end

%startyf = XYTsolute_cell{1,2};
startxf = zeros(size(startyf));

widths = zeros(size(XYTsolute_cell{1,3},2)-1,1);
xMaxes = zeros(size(XYTsolute_cell{1,3},2)-1,1);
xLefts = zeros(size(XYTsolute_cell{1,3},2)-1,1);
xRights = zeros(size(XYTsolute_cell{1,3},2)-1,1);
maxms = zeros(size(XYTsolute_cell{1,3},2)-1,1);

flux_inlet = zeros(size(XYTsolute_cell{1,3},2)-1,1);

time = XYTcolloid_cell{1,3};
num_particles = zeros(1,size(time,2));

count = 1;
while count <= size(time,2)

% march in time
% velocity evaluate at cell boundaries, concentration in center
[vp_x,vp_y,vf_x,vf_y,dcdx,dcdy,Vs,dVsdx]=vp_discrete(c_interpolated(:,:,count),...
    XYTsolute_cell{1,1},XYTsolute_cell{1,2},colloid_log.dx,...
    colloid_log.dy,colloid_log.Gp_Ds,solute_log.Gw_Ds,solute_log.L_h);
    
[dndx,dndy]=n_flux(colloid_cell{count,1},...
    XYTcolloid_cell{1,1},XYTcolloid_cell{1,2},colloid_log.dx,...
    colloid_log.dy);
% the X=0 values destroy the scale
%dndx(1,:)=0; dndy(1,:)=0;

n = colloid_cell{count,1};
y = XYTcolloid_cell{1,2};
y_c=0.5*(y(2:end)+y(1:end-1));
f_l = zeros(1,length(y_c));
j=1:length(y_c);
vp_l = (vp_x(1,j)+vp_x(1,j+1))*0.5;
n_l = n(1,j);
f_l(1,j) = vp_l.*n_l - colloid_log.Dp_Ds*(n(2,j)-n(1,j))/colloid_log.dx;
flux_inlet(count) = trapz(y_c, f_l(1,:));


xxx = XYTcolloid_cell{1,1}'; nnn = mean(colloid_cell{count,1},2);
% cut off first index bc instability at x=0
[widths(count),xMaxes(count),xLefts(count),xRights(count),maxms(count)]...
    = width_nAvg(xxx(2:end),nnn(2:end));

int_x = trapz(XYTcolloid_cell{1,1},colloid_cell{count,1},1);
int_xy = trapz(XYTcolloid_cell{1,2},int_x);
num_particles(count) = int_xy; %#ok<SAGROW>


count = count + 1;
end


cumulative_nParticles = trapz(time(1:time_integrate_range)...
    ,num_particles(1:time_integrate_range));


[xx_cCell,yy_cCell]=meshgrid(XYTsolute_cell{1,1},XYTsolute_cell{1,2});
[xx_nCell,yy_nCell]=meshgrid(XYTcolloid_cell{1,1},XYTcolloid_cell{1,2});

function [dndx,dndy]=n_flux(n,x,y,dx,dy)

dndx = [(diff(n)/dx); zeros(1,length(y))];
dndy = [(diff(n')/dy)', zeros(length(x),1)];

end


