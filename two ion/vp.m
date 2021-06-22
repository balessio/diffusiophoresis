function [vp_x,vp_y,vf_x,vf_y,dcdx,dcdy]=vp(x,y,t,GpDs,GwDs,B,L_h)

% x is a row vector
% y is a column vector 

[xx,yy]=meshgrid(x,y);

res=1.0;
c=zeros(size(xx));
i=0;
dcdx=zeros(size(xx)); dcdy=zeros(size(xx));
dc2dx2=zeros(size(xx));

while(res>1.0e-4)
lambda=(2*i+1)*pi/2;
val=c;
c = c + (1-B)*2/lambda*sin(lambda*xx)*exp(-lambda^2*t);
res=sqrt(sum(sum((c-val).^2))./sum(sum(c.^2)));
dcdx = dcdx + (1-B)*2*cos(lambda*xx)*exp(-lambda^2*t);
dc2dx2 = dc2dx2 -(1-B)*2*lambda*sin(lambda*xx)*exp(-lambda^2*t);
i=i+1;
end
c=c+B;

vf_x = -1/2*GwDs*dcdx./c.*(3*yy.^2*L_h^2-1);
vf_y = 1/2*GwDs*yy.*(yy.^2*L_h^2-1).*(dc2dx2./c - (dcdx./c).^2);

vDp_x = GpDs*dcdx./c;
vDp_y = 0;

vp_x= vf_x + vDp_x;
vp_y = vf_y + vDp_y;
end
