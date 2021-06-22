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
dt = parameter_file_unload.dt;
output_folder = parameter_file_unload.output_folder;

Gw_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiW+PsiW^2/8);
Gp_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiD+PsiD^2/8);

Dp_Ds = Dp_D_Na;

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

save(join([output_folder,'log.mat'],''),'Gp_Ds','Gw_Ds','Dp_Ds','L_h','B',...
'n_bulk','dx','dt','t0','duration');

% location of cell boundaries
x=0:dx:1;
t=t0:dt:duration;

% location of cell centers
x_c=0.5*(x(2:end)+x(1:end-1));

% initialize n
n = exp(-(x_c-0.2).^2*100);

f_l=n.*0;
f_r=n.*0;

XYTcolloid_cell=cell(1,2);
XYTcolloid_cell{1,1}=x_c; XYTcolloid_cell{1,2}=t;
colloid_cell=cell(length(t),1);

tic;
% march in time
for count=1:length(t)
    
    % velocity evaluate at cell boundaries, concentration in center
    %[V_DP,V_s] = vp_dispersion(x,t(count),Gp_Ds,Gw_Ds,B);
    
    n0=n;
    
    n=fsolve(@(Xval)myfun(Xval,x,x_c,dx,t(count),dt,B,n_bulk,...
        Gp_Ds,Gw_Ds,Dp_D_Na,L_h,n0),reshape(n0,[length(x_c),1]));
    n=reshape(n, [length(x_c),1]);
        
    colloid_cell{count,1} = n;
    
    figure(1)
    drawnow
    
    plot(x_c,n,'-k','linewidth',1.5)
    %axis([0 1 0 max(n)*1.25])
    axis([0 1 0 1.0])
    %title('Y-averaged colloid concentration vs X')
    title('y-averaged colloid concentration vs X, tau='+string(t(count)))
    xlabel('X')
    ylabel('< n >')
    grid on

end
toc;

save(join([output_folder,'XYTcolloid_cell.mat'],''), 'XYTcolloid_cell');
save(join([output_folder,'colloid_cell.mat'],''), 'colloid_cell');

function [F,f_l,f_r]=myfun(Xval,x,x_c,dx,t,dt,B,n_bulk,GpDs,GwDs,DpDs,L_h,n0)

n=reshape(Xval,[length(x_c),1]);

[V_DP,V_s] = vp_dispersion(x_c,t,GpDs,GwDs,B);
V_dp = transpose(V_DP); V_S = transpose(V_s);

% left flux
i=2:length(x_c);
n_l = (V_dp(i)>0).*n(i-1) + (V_dp(i) <= 0).*n(i);
f_l(i) = V_dp(i).*n_l - ((DpDs + (2/105)*L_h^(-2)*DpDs^(-1)*V_S(i).^2)).*(n(i)-n(i-1))/dx;
%f_l(1)= - ((DpDs + (2/105)*L_h^(-2)*DpDs^(-1)*V_S(1)^2))*(n(1)-n_bulk)/dx/2;
f_l(1)= - ((DpDs ))*(n(1)-n_bulk)/dx/2;

% right flux
i=1:length(x_c)-1;
n_r = (V_dp(i+1)>0).*n(i) + (V_dp(i+1) <= 0).*n(i+1);
f_r(i) = V_dp(i+1).*n_r - ((DpDs + (2/105)*L_h^(-2)*DpDs^(-1)*V_S(i+1).^2)).*(n(i+1)-n(i))/dx;
f_r(length(x_c))=0;

F = (n-n0 - dt*((f_l-f_r)/dx ));

end
