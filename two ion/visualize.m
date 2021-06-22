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

count = ind_closest(which_time,XYTcolloid_cell{1,3});

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

num_particles = zeros(1,length(XYTcolloid_cell{1,3}));

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
    
[xx_cCell,yy_cCell]=meshgrid(XYTsolute_cell{1,1},XYTsolute_cell{1,2});
[xx_nCell,yy_nCell]=meshgrid(XYTcolloid_cell{1,1},XYTcolloid_cell{1,2});
    

function [dndx,dndy]=n_flux(n,x,y,dx,dy)

dndx = [(diff(n)/dx); zeros(1,length(y))];
dndy = [(diff(n')/dy)', zeros(length(x),1)];

end
