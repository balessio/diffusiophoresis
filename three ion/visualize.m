%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global colloid_folder
global which_time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colloid_cell_unload = matfile(join([colloid_folder,'colloid_cell.mat'],''));
colloid_cell = colloid_cell_unload.colloid_cell;
XYTcolloid_cell_unload = matfile(join([colloid_folder,'XYTcolloid_cell.mat'],''));
XYTcolloid_cell = XYTcolloid_cell_unload.XYTcolloid_cell;
time_n = XYTcolloid_cell{1,3}';

potentials_unload = matfile(join([colloid_folder,'potentials.mat'],''));
PsiW_array = potentials_unload.PsiW_array;
PsiD_array = potentials_unload.PsiD_array;

colloid_log = matfile(join([colloid_folder,'log.mat'],''));
solute_log = matfile(join([colloid_log.solute_cell_folder,'log.mat'],''));
solute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'solute_cell.mat'],''));
solute_cell = solute_cell_unload.solute_cell;
XYTsolute_cell_unload = matfile(join([colloid_log.solute_cell_folder,'XYTsolute_cell.mat'],''));
XYTsolute_cell = XYTsolute_cell_unload.XYTsolute_cell;
time_c = XYTsolute_cell{1,3}';

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

count = ind_closest(which_time,time_c);

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dndx,dndy] = n_flux(n,x,y,dx,dy)
dndx = [(diff(n)/dx); zeros(1,length(y))];
dndy = [(diff(n')/dy)', zeros(length(x),1)];
end

