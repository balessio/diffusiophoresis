function [vp_x,vp_y,vf_x,vf_y,vDp_x,vDp_y]=vp_discrete_multi_ion(...
    c1,c2,c3,Di,zi,x,y,dx,dy,Pe,PsiW,PsiD,L_h)

% x is a row vector
% y is a column vector
[xx,yy]=meshgrid(x,y);

dc1dx = [(c1(2,:)-c1(1,:))/(2*dx);...
    (c1(3:end,:)-c1(1:end-2,:))/(2*dx);zeros(1,length(y))];
dc1dy = [zeros(length(x),1),...
    (c1(:,3:end)-c1(:,1:end-2))/(2*dy),zeros(length(x),1)];
dc2dx = [(c2(2,:)-c2(1,:))/(2*dx);...
    (c2(3:end,:)-c2(1:end-2,:))/(2*dx);zeros(1,length(y))];
dc2dy = [zeros(length(x),1),...
    (c2(:,3:end)-c2(:,1:end-2))/(2*dy),zeros(length(x),1)];
dc3dx = [(c3(2,:)-c3(1,:))/(2*dx);...
    (c3(3:end,:)-c3(1:end-2,:))/(2*dx);zeros(1,length(y))];
dc3dy = [zeros(length(x),1),...
    (c3(:,3:end)-c3(:,1:end-2))/(2*dy),zeros(length(x),1)];

vs_n1 = Di(1)*zi(1)*dc1dx(:,1) + Di(2)*zi(2)*dc2dx(:,1) + Di(3)*zi(3)*dc3dx(:,1);
vs_d1 = Di(1)*zi(1)^2*c1(:,1) + Di(2)*zi(2)^2*c2(:,1) + Di(3)*zi(3)^2*c3(:,1);
vs_n2 = zi(1)^2*dc1dx(:,1) + zi(2)^2*dc2dx(:,1) + zi(3)^2*dc3dx(:,1);
vs_d2 = zi(1)^2*c1(:,1) + zi(2)^2*c2(:,1) + zi(3)^2*c3(:,1);
Vs_calc = PsiW(:,1).*vs_n1./vs_d1 + (PsiW(:,1).^2/8).*vs_n2./vs_d2;
dVsdx_calc = [(Vs_calc(2)-Vs_calc(1))/(2*dx);...
    (Vs_calc(3:end)-Vs_calc(1:end-2))/(2*dx);0];
Vs = zeros(size(xx));
dVsdx = zeros(size(xx));
for i=1:length(y)
    Vs(i,:)=Pe*Vs_calc';
    dVsdx(i,:)=Pe*dVsdx_calc';
end

vf_x = -1/2*Vs.*(3*yy.^2*L_h^2-1);
vf_y = 1/2*yy.*(yy.^2*L_h^2-1).*dVsdx;

vdpx_n1 = Di(1)*zi(1)*dc1dx + Di(2)*zi(2)*dc2dx + Di(3)*zi(3)*dc3dx;
vdpx_d1 = Di(1)*zi(1)^2*c1 + Di(2)*zi(2)^2*c2 + Di(3)*zi(3)^2*c3;
vdpx_n2 = zi(1)^2*dc1dx + zi(2)^2*dc2dx + zi(3)^2*dc3dx;
vdpx_d2 = zi(1)^2*c1 + zi(2)^2*c2 + zi(3)^2*c3;
vDp_x = PsiD.*vdpx_n1./vdpx_d1 + (PsiD.^2/8).*vdpx_n2./vdpx_d2;

vdpy_n1 = Di(1)*zi(1)*dc1dy + Di(2)*zi(2)*dc2dy + Di(3)*zi(3)*dc3dy;
vdpy_d1 = Di(1)*zi(1)^2*c1 + Di(2)*zi(2)^2*c2 + Di(3)*zi(3)^2*c3;
vdpy_n2 = zi(1)^2*dc1dy + zi(2)^2*dc2dy + zi(3)^2*dc3dy;
vdpy_d2 = zi(1)^2*c1 + zi(2)^2*c2 + zi(3)^2*c3;
vDp_y = PsiD.*vdpy_n1./vdpy_d1 + (PsiD.^2/8).*vdpy_n2./vdpy_d2;

vDp_x = Pe*vDp_x'; vDp_y = Pe*vDp_y';

vp_x= vf_x + vDp_x;
vp_y = vf_y + vDp_y;

end
