function [vp_x,vp_y,vf_x,vf_y,dcdx,dcdy,Vs,dVsdx]=vp_discrete(c,x,y,dx,dy,Gp_Ds,Gw_Ds,L_h)

% x is a row vector
% y is a column vector
[xx,yy]=meshgrid(x,y);

%dcdx = [(diff(c)/dx); zeros(1,length(y))];
%dcdy = [(diff(c')/dy)', zeros(length(x),1)];

dcdx = [(c(2,:)-c(1,:))/(2*dx);...
    (c(3:end,:)-c(1:end-2,:))/(2*dx);zeros(1,length(y))];
dcdy = [zeros(length(x),1),...
    (c(:,3:end)-c(:,1:end-2))/(2*dy),zeros(length(x),1)];

c_wall = c(:,1);
dc_walldx = dcdx(:,1);
%d2c_walldx2 = [(diff(dc_walldx)/dx); 0];
d2c_walldx2 = [(c(3,1)-2*c(2,1)+c(1,1))/dx^2;...
    (c(3:end,1)-2*c(2:end-1,1)+c(1:end-2,1))/dx^2; 0];
%exp_decay = 100;
%Vs_calc = -Gw_Ds*dc_walldx.*(1-exp(-exp_decay*x))'./c_wall;
%dVsdx_calc = -Gw_Ds*(d2c_walldx2./c_wall - (dc_walldx./c_wall).^2).*(1-exp(-exp_decay*x))';
Vs_calc = -Gw_Ds*dc_walldx./c_wall;
dVsdx_calc = -Gw_Ds*(d2c_walldx2./c_wall - (dc_walldx./c_wall).^2);
Vs = zeros(size(xx));
dVsdx = zeros(size(xx));
for i=1:length(y)
    Vs(i,:)=Vs_calc';
    dVsdx(i,:)=dVsdx_calc';
end
vf_x = 1/2*Vs.*(3*yy.^2*L_h^2-1);
vf_y = -1/2*yy.*(yy.^2*L_h^2-1).*dVsdx;

Vs = Vs(1,:); dVsdx = dVsdx(1,:);

vDp_x = Gp_Ds*dcdx'./c';
vDp_y = Gp_Ds*dcdy'./c';

vp_x= vf_x + vDp_x;
vp_y = vf_y + vDp_y;
%vp_y = vDp_y;

end
