function [V_DP,V_s]=vp(x,t,GpDs,GwDs,B)

% x is a row vector

res=1.0;
c=zeros(size(x));
i=0;
dcdx=zeros(size(x)); dc2dx2=zeros(size(x));

while(abs(res)>1.0e-4)
lambda=(2*i+1)*pi/2;
val=c;
c = c + (1-B)*2/lambda*sin(lambda*x)*exp(-lambda^2*t);
res=sqrt(sum(sum((c-val).^2))./sum(sum(c.^2)));
dcdx = dcdx + (1-B)*2*cos(lambda*x)*exp(-lambda^2*t);
dc2dx2 = dc2dx2 -(1-B)*2*lambda*sin(lambda*x)*exp(-lambda^2*t);
i=i+1;
end
c=c+B;

V_DP = GpDs*dcdx./c;
V_s = GwDs*dcdx./c;


end
